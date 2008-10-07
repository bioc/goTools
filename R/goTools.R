######################
## New version of goTools function compatible with
## new version of GO annotation package (GO.db)
## Author: Agnes Paquet
## Date: October 4, 2008


##source("C:/MyDoc/Projects/madman/Rpacks-devel/goTools/R/goTools_V3.R")
## library(GO.db)
## 

ontoCompare <- function(genelist,probeType=c("GO","hgu133a"),goType="All",
                        endnode, method=c("TGenes","TIDS","none"),plot=FALSE,...)
  {
    cat("Starting ontoCompare...\n")
    probeType <- probeType[1]
    if(missing(endnode))
      endnode <- EndNodeList()
    else
       endnode <- unique(c("GO:0003673","GO:0003674", "GO:0005575", "GO:0008150",
                           endnode))

    if(missing(method))
      method <- "TGenes"
    else
      method <- method[1]

    FullGOenv <- as.list(GOTERM)
    if (probeType!="GO")
      {
        ## Convert the probe ids to GO
        golist <- lapply(genelist,getGOID,probeType=probeType)
      }

    else
      {
        golist <- genelist
      }

    golistOK <- c()
    for(i in 1:length(golist))
          {
            tmp <- golist[[i]]
            tmp2 <- lapply(tmp,function(x){return(x[x %in% names(FullGOenv)])})
            golistOK <- c(golistOK,list(tmp2))
          }
    
    names(golistOK) <- names(golist)
    
    ##Now we have a list of list of goIds
    ## Get the corresponding endnodes

    parentsIDs <- c()
    for(i in 1:length(golistOK))
      {
        tmp <- golistOK[[i]]
        tmp2 <- lapply(tmp,parentsVectWraper,endnode)
        parentsIDs <- c(parentsIDs, list(tmp2))
      }
    names(parentsIDs) <- names(golistOK)

 
    finalParentsIDs <- c()
    ## Keep only Ontology of interest as specified in goType
    for(i in 1:length(parentsIDs))
      {
        tmp <- parentsIDs[[i]]
        tmp2 <- lapply(tmp,getOntology,goType,FullGOenv)
        finalParentsIDs <- c(finalParentsIDs,list(tmp2))
      }

    ## Reference endnode for the plot
    refnode <- getOntology(endnode,goType,FullGOenv)
    
    ## Now that we have the mapping for the parents
    ## Count how many

    results <- ontoCompare.main(finalParentsIDs,method,refnode,FullGOenv)
    if(plot)
      ontoPlot(results,...)
    return(results)
  }


ontoCompare.main <- function(parentslist,method=c("TGenes","TIDS","none"),refnode,FullGOenv=NULL)
  {
    method <- method[1]
    if(missing(FullGOenv))
      FullGOenv <- as.list(GOTERM)
    resmat <- c()
    for(i in 1:length(parentslist))
      {
        tmp <- parentslist[[i]]
        nb <- table(unlist(tmp))
        index <- match(names(nb),refnode)
        nas <- sum(is.na(unlist(tmp)))
        nb <- c(nb,nas)
        res <- rep(0,length(refnode)+1)
        if(method=="TGenes")
          {
            den <- length(tmp)
            res[c(index,length(res))] <- nb/den
          }
        if(method=="TIDS")
          {
            den <- length(unlist(tmp))
            res[c(index,length(res))] <- nb/den
          }
        if(method=="none")
          {
            res[c(index,length(res))] <- nb
          }
        resmat <- cbind(resmat,res)
      }
    labs <- sapply(FullGOenv[refnode],Term)
    rownames(resmat) <- c(labs,"NotFound")
    colnames(resmat) <- names(parentslist)
    ##return(resmat)
    return( as.matrix(resmat[apply(resmat,1,function(x){sum(x)>0}),]))
  }



getOntology <- function(idvect,goType=c("All","BP","CC","MF"), FullGOenv)
  {
    if(is.na(idvect[1]))
      {
        return(NA)
      }
    else
      {
        goType=goType[1]
        if(goType=="All")
          return(idvect)
        else
          {
            index <- match(idvect,names(FullGOenv))
            cat <-  sapply(FullGOenv[index],Ontology)
            return(idvect[cat==goType])
          }
      }
  }

EndNodeList <- function() {
  MFendnode <- get("GO:0003674", env=GOMFCHILDREN)
  CCendnode <- get("GO:0005575", env=GOCCCHILDREN)
  BPendnode <- get("GO:0008150", env=GOBPCHILDREN)
  EndNodeList <- c("GO:0003674", "GO:0005575", "GO:0008150",
                   MFendnode, CCendnode, BPendnode)
  return(EndNodeList)
}

## Take probe ids, return GO, for any available platform
getGOID <- function (x, probeType="operon")
{
  if(probeType == "operon")
    print("operon is not supperted")
    ## nothing for now
    ## res <- getGO.operon(x)
  else
    {
      library(paste(probeType,".db",sep=""), character.only = TRUE)
      GOenv <- get(paste(probeType, "GO", sep = ""))
      tmp <- mget(x, env = GOenv)
      res <- lapply(tmp, names)
    }
  return(res)
}

isEndNode <- function(id, endnode) {
  if(missing(endnode))
    endnode <- EndNodeList()
  return(is.element(id, endnode))

}

parentsVectWraper <- function(goidvect,endnode)
  {
    if (is.null(goidvect))
      {
        results <- NA
      }
    else
      {
        results <- c()
        test <- sapply(goidvect, isEndNode,endnode=endnode)
        results <- c(results,goidvect[test]) ## Keep a list of the endnode we hit
        var <- goidvect[!test]
        while(length(var) !=0)
          {
            parents <- unlist(sapply(var,goParents))
            test <- unlist(lapply(parents,isEndNode, endnode))
            results <- c(results,parents[test])
            var <- parents[!test]
          }
        if(length(results)==0)
          results <- NA
      }
    return(unique(results))
  }

goParents <- function(goid)
  {
    cat <- Ontology(get(goid,GOTERM))
    if(cat == "all")
      {
        stop("reach the top of the tree")
      }
    else
      {
        envi <- get(paste("GO",cat,"PARENTS",sep=""),
                    env=as.environment("package:GO.db"))
        res <- get(goid,env=envi)
      }
    return(res)
  }

#####################################################
# Function:ontoPlot
#####################################################
## Suppose we have a list  resultsing from getGOList

ontoPlot <- function(objM,
                     names.arg=NULL,
                     beside=TRUE,
                     las=2,
                     legend.text=TRUE,
                     ...)
  {
    ## obj <- results from ontoCompare
    names.arg=rownames(objM)
    if(ncol(objM) == 1)
      pie(as.vector(objM[objM!=0]), labels=names.arg, ...)
    else {
      if(is.logical(legend.text))
        {
          if(legend.text)
            legend.text <- colnames(objM)
          else
            legend.text <- NULL
        }
      x <- barplot(t(objM), col=rainbow(ncol(objM)), beside=beside, las=las,
                   names.arg=names.arg,legend.text=legend.text, ...)
      return(x)
    }
  }

###########################################################
## Other functions from previous veriosn of goTools
## May not work...

goChildren <- function(id) {
  ## Assume id is valid.
  cat <- Ontology(get(id,GOTERM))
  if(!is.na(id) & !setequal(cat, "GO"))
    {
      envi <- get(paste("GO", cat, "CHILDREN", sep=""),
                  envir=as.environment("package:GO.db"))
      if (id %in% names(as.list(envi)))
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
      cat(paste("rank=",i,"\n"))
      cust <- unique(unlist(lapply(cust, goChildren)))
      cust <- cust[!is.na(cust)]
      res <- c(res,cust)
    }
  return(unique(res))
}
