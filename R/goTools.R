#####################################################
# Gotools.R
#
# Modified: April 11, 2005
# Functions for description of oligos using Gene Ontology
#
# TO USE
#
# ontoCompare(list(id, id), method="TIDS", probeType="rgu34a", goType="BP", plot=TRUE)
#
#####################################################


#####################################################
# Function: getGOCategory
#####################################################
# This function returns the GO category corresponding
# to one id. Needs to be entered as "GO:0000000"


getGoCategory <- function(id) {
  cat <- Ontology(get(id, env=GOTERM))
  return(cat)
}


#####################################################
# Function: getGOTerm
#####################################################

getGOTerm <- function (num)
{
  if (!nchar(num[1]))
    return(list())
  if (exists(num, GOTERM))
    {
      res <- get(num, env=GOTERM)

      if (Ontology(res)== "BP")
        return(list(name = Term(res), type = "Biological Process"))
      if (Ontology(res)== "CC")
        return(list(name = Term(res), type = "Cellular Component"))
      if(Ontology(res) == "MF")
        return(list(name = Term(res), type = "Molecular Function"))
    }
}




#####################################################
# Function: getGOList
#####################################################

getGOList <- function(numvect, goType=c("All", "BP", "CC", "MF"))
{
  ## print("in getGOList")
  ## numvect is a vector
  results <- NA
  if(!is.null(numvect))
    {
      if(!is.na(numvect[1]))
        {
          res <- sapply(numvect, getGOTerm)
          if(length(res) !=0)
            {
              goType <- goType[1]
              if(goType=="All")
                results <- res
              if(goType=="BP")
                results <- as.matrix(res[, unlist(res["type",]) == "Biological Process"])
              if(goType=="CC")
                results <- as.matrix(res[, unlist(res["type",]) == "Cellular Component"])
              if(goType=="MF")
                results <- as.matrix(res[, unlist(res["type",]) == "Molecular Function"])
              if(ncol(results) == 0)
                results <- NA
            }}}
  return(results)
}

#####################################################
# Function goParents
#####################################################
#Given a GO id "GO:XXXXXX", returns its parents

goParents <- function(id) {
  ## Assume id is valid.
  cat <- getGoCategory(id)
  if (cat == "GO")
    stop("We have reached the top node of the GO tree")
  else
    {
      envi <- get(paste("GO", cat, "PARENTS", sep=""),
                  envir=as.environment("package:GO"))
      res <- get(id, env=envi)
    }
  return(res)
}




#####################################################
# Function EndNodeList
#####################################################
# Returns the GO end-nodes before MF, BP, CC
# as a look-up vector

EndNodeList <- function() {
MFendnode <- get("GO:0003674", env=GOMFCHILDREN)
  CCendnode <- get("GO:0005575", env=GOCCCHILDREN)
  BPendnode <- get("GO:0008150", env=GOBPCHILDREN)
  EndNodeList <- c("GO:0003674", "GO:0005575", "GO:0008150",
                   MFendnode, CCendnode, BPendnode)
  return(EndNodeList)
}

#####################################################
# Function goChildren
#####################################################
# Returns the children of a GO id
# BUG FIX Jean: don't need to check for the top node for this function

goChildren <- function(id) {
## Assume id is valid.
  cat <- getGoCategory(id)
  if(!is.na(id) & !setequal(cat, "GO"))
    {
      envi <- get(paste("GO", cat, "CHILDREN", sep=""),
                  envir=as.environment("package:GO"))
      if (id %in% ls(env=envi))
        {
          res <- get(id, env=envi)
          return(res[!is.na(res)])
        }
      else
        return(NA)
    }
  else
    return(NA)
}


#####################################################
# Function CustomEndNodeList
#####################################################
# Create a list of all end-nodes going down the GO tree
# from the top to rank

CustomEndNodeList <- function(id,rank=1){
  cust <- id
  res <- numeric(0)
  for(i in 1:rank)
    {
      print(paste("rank=",i))
      cust <- unique(unlist(lapply(cust, goChildren)))
      cust <- cust[!is.na(cust)]
      res <- c(res,cust)
    }
  return(unique(res))
}


#####################################################
# Function: parentsList
#####################################################
#  Returns a list of parents, given a list of Go ids
#  Input is a list or vector, output a vect
## Assume all GO are in the metadata
parentsList <- function(vect) {
  
  if(is.list(vect))
    pars <- lapply(vect, sapply, goParents)
  else
    pars <- sapply(vect, goParents)
    
  pars <- unlist(pars)
  return(as.vector(pars))
}

#####################################################
# Function: isEndNode
#####################################################
#Test if a given id is in the look-up table

isEndNode <- function(id, endnode) {
  if(missing(endnode))
    endnode <- EndNodeList()
  return(is.element(id, endnode))

}

#####################################################
# Function: goId
#####################################################

getGOID <- function (x, probeType="operon")
{
  #print("in getGOID")
  if(probeType == "operon")
    res <- getGO.operon(x)
  else
    {
      library(probeType, character.only = TRUE)
      GOenv <- get(paste(probeType, "GO", sep = ""))
      tmp <- mget(x, env = GOenv)
      res <- lapply(tmp, names)
    }
  return(res)
}

getGO.operon <- function(oligo, gotableinput)
  {
    if(missing(gotableinput))
      if(!("gotable" %in% ls(1)))
        assign("gotable", updateOligo2GO(), envir=.GlobalEnv)

    if(!missing(gotableinput))
      gotable <- gotableinput
    #print("here")
    print(dim(gotable))#not doing anything?  Mainly for debug
    res <- lapply(oligo, getGO.operon.main, gotable=gotable)
    names(res) <- oligo

    index <- sapply(res, function(z) {
      if(length(z)==1)
        {
          if(is.na(z)) return(FALSE)
          else return(TRUE)
        }
      else return(TRUE)
    })
    res <- res[index]
                  
    return(res)
  }

getGO.operon.main <- function(oligo, gotable)
  {
    vect <- NA
    #print("there")
    ind <- grep(oligo, as.character(gotable[,1]))
    #print("there 2")
    if (!setequal(ind, numeric(0))) {
      vect <- unlist(strsplit(as.character(gotable[ind,2]), split=" :: "))
      if(!is.null(vect)) vect <- gsub(" ", "", vect)
      #vv <- vect[vect %in% ls(GOCATEGORY)]
      vv <- vect[vect %in% ls(env = GOTERM)]
      vv2 <- vv[!is.na(sapply(vv,getGoCategory))]
      vect <- vv2
    }
    return(vect)
  }

#####################################################
# Function:goBarBarplot
#####################################################
## Suppose we have a list  resultsing from getGOList

ontoPlot <- function(objM,
                        beside=TRUE,
                        las=2,
                        legend.text=TRUE,
                        ...)
  {
    ## obj <- results from ontoCompare
    if(nrow(objM) == 1)
      pie(as.vector(objM[objM!=0]), labels=colnames(objM), ...)
    else {
      if(is.logical(legend.text))
        {
          if(legend.text)
            legend.text <- rownames(objM)
          else
            legend.text <- NULL
        }
      x <- barplot(objM, col=rainbow(nrow(objM)), beside=beside, las=las,
                   legend.text=legend.text, ...)
      return(x)
    }
  }


#####################################################
# Function: updateOligo2GO
#####################################################
## Function to pull annotations from the web site

## Still reloading the table each time,
## unless you save it to the env before
updateOligo2GO <- function(url)
  {
    cat("Downloading Oligo 2 GO annotation table ...")
    if(missing(url))
      url <- "http://arrays.ucsf.edu/download/GO-IDs"
    gotable <- read.table(url, header=TRUE, sep="\t", fill=TRUE)
    ind <- match(unique(gotable[,1]), gotable[,1])
    gotable <- gotable[ind,]
    cat("done.\n")
    return(gotable)
  }

#################################################
## Function gowraper
#################################################
## Given oligo id, returns its end-nodes

gowraper <- function(oligo, endnode, probeType)
  {
   # print("in gowraper")
    if(missing(endnode))
      endnode <- EndNodeList()

    if(missing(probeType))
      probeType <- "operon"

    #print(paste("probType= ", probeType))
    goItmp <- getGOID(oligo, probeType=probeType)
    ## goItmp is a list of GO for each oligo ID

    ## Check go exists in data base ## It should but version differences
    #FULLGOList <- ls(GOCATEGORY)
    FULLGOList <- ls(env = GOTERM)  ## List of all GOTERM

    goItmp2 <- lapply(goItmp, function(x){x[x %in% FULLGOList]})
    goI <- lapply(goItmp2, function(x){x[!is.na(sapply(x,getGoCategory))]})

    #goI <- lapply(goItmp, function(x){x[x %in% FULLGOList & !is.na(sapply(x,getGoCategory))]})
    ## remove all the names that are not in GOTERM
    
    ## Find parents
    results <- lapply(goI, parentsVectWraper, endnode)
    ## List of vector: names(results) = oligo ID and each
    ## vector represent the end node GO term. 
    
    #print("end gowraper")
    return(results)
  }

#################################################
## Function parentsListWraper
#################################################
## Given a list of GO ids and a list of GO endnodes,
## returns the GO endnodes parents of the GO ids

parentsListWraper <- function(goI, endnode, listres = TRUE)
  {
    ## input goI is a list
    if(missing(endnode))
      endnode <- EndNodeList()
    results <- lapply(goI, parentsVectWraper, endnode)
    if(listres)
      return(results)
    else
      return(unique(unlist(results)))
  }

parentsVectWraper <- function(goI, endnode)
  {
    #print("in parentsVectWraper")
    ## input goI is a vect
    if(missing(endnode))
      endnode <- EndNodeList()

    if(is.null(goI))
      {
        results <- NA
      }
    else
      {
        if(!is.na(goI[1]))
          {
            results <- NULL
            test <- sapply(goI, isEndNode, endnode=endnode)
            results <-  c(results, goI[test])
            var <- goI[!test]
            while (length(var) != 0)
               {
                 parents <- parentsList(var)
                 test <- unlist(lapply(parents, isEndNode, endnode))
                 results <- c(results, parents[test])
                 var <- parents[!test]
               }
            if(is.null(results)) results <- NA
          }
        else
          {
            results <- NA
          }
      }
    return(unique(results))
  }


#################################################
## Function
#################################################
## Given oligo ids, return percentage of each
## end-node, and the % of gene without annotations
## data describe our object.  It can be
## probeType == "operon": List of operon ID
## probeType == "hgu133a" : List of Affy ID from chip hgu133a
## probeType == "GO" : we have GO ids e.g. GO:xxxxxx
## probeType == "GONames" : We have GO description e.g. "transcription factor

ontoCompare  <- function(obj,  method=c("TGenes", "TIDS", "none"),
                         probeType=c("GO", "operon"), goType="All", plot=FALSE,
                         endnode, ...)
  {

    print("Starting ontoCompare...")
    probeType <- probeType[1]
    if(missing(endnode))
      endnode <- EndNodeList()

    ## We need to add GO, Mf, BP, CC in the list, to stop the recurrence
    endnode <- unique(c("GO:0003673","GO:0003674", "GO:0005575", "GO:0008150",endnode))

    ## List of GO or list of oligo

    if(probeType == "GO")
      {
        GOID <- list()
        for(i in obj)
          {
            goItmp <- i
            FULLGOList <- ls(env = GOTERM)  ## List of all GOTERM
            goItmp2 <- lapply(goItmp, function(x){x[x %in% FULLGOList]})
            goI <- lapply(goItmp2, function(x){x[!is.na(sapply(x,getGoCategory))]})

            #goI <- lapply(goItmp,
            #              function(x){x[x %in% FULLGOList & !is.na(sapply(x,getGoCategory))]})
            ## remove all the names that are not in GOTERM
            ## Find parents
            GOID <- c(GOID, list(lapply(goI, parentsVectWraper, endnode)))
          }
        names(GOID) <- names(obj)
      }
    else
      {
        GOID <- lapply(obj, gowraper, endnode=endnode, probeType=probeType)
    }
    ## GOID = GO 
    objlist <- lapply(GOID, function(x){lapply(x, getGOList, goType=goType)})
   
    ## List of GONames: Description

    res <- ontoCompare.main(objlist, method = method[1])

    if(plot)
      ontoPlot(res, ...)

    return(res)
  }

## Given oligo ids, return percentage of each
## end-node, and the % of gene without annotations
ontoCompare.main <- function(obj, method=c("TGenes", "TIDS", "none"))
  {
    ## obj is a list of list.  So if you only have one element,
    ## you still need to create a list of length 1.
    method <- method[1]
    if(length(obj) > 1)
      {
        #print("ontoCompare 3")
        NotFoundGenes <-  unlist(lapply(obj, function(x){sum(is.na(unlist(x)))}))
        newobj <- lapply(obj, function(x){
          y <- lapply(x[!is.na(x)], function(x){unlist(x[1,])})
          table(unlist(y))})

        print("Number of lists > 1")
        x <- unique(unlist(lapply(newobj, function(x){return(names(x))})))
        newM <- matrix(0, ncol=length(x) + 1, nrow=length(newobj))
        colnames(newM) <- c(x, "NotFound")
        newM[,length(x) + 1] <- NotFoundGenes

        if(is.null(names(obj)))
          rownames(newM) <- as.character(1:length(obj))
        else
          rownames(newM) <- names(obj)

        ## Enter probeType
        for(i in 1:length(newobj))
          {
            newM[i, names(newobj[[i]])] <- newobj[[i]]
          }
      }
    else
      {
        print("Number of lists = 1")
        obj2 <- obj[[1]]
        xx <- table(unlist(lapply(obj2[!is.na(obj2)], function(x){unlist(x[1,])})))
        newM <- matrix(c(xx, sum(is.na(unlist(obj2)))), nrow=1)
        rownames(newM) <- "1"
        colnames(newM) <- c(names(xx), "NotFound")
      }

    print(paste("Using method:",method[1]))
    if(method == "TGenes")
      {
        TGenes <- unlist(lapply(obj, length))
        res <- newM / TGenes
      }
    if(method=="TIDS")
      {
        TIDs <- lapply(lapply(obj, unlist), length)
        res <- newM / unlist(TIDs)
      }
    if(method=="none")
      res <- newM
    return(res)
  }

