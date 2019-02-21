#
# GET CO-EXPRESSION NETWORK 
getTopCoexprPairsBy <- function(exprMat, cutOff, self = F, mode="quantile", debug=F) {
  #remove zero expressions
  if (debug) print("removing genes expressing less than one FPKM or TPM")
  
  minValue=-2
  if (debug) print("correlation generating")
  if (mode=="quantile") {
    corMat = cor(t(exprMat))
    #if(debug) print(corMat)
    
    if(debug) print("na amount:")
    if(debug) print(table(is.na(corMat)))
    
    corMat[is.na(corMat)]=minValue
    cutOff = as.numeric( quantile(corMat[lower.tri(corMat)], cutOff, na.rm=T) )
    
    if(debug) print("quantile: ")
    if(debug) print(cutOff)
    
    resTable = melt(corMat)[,1:3]  
    #resTable[,1] = sapply(resTable[,1], as.character)
    #resTable[,2] = sapply(resTable[,2], as.character)
    if(debug) print(class(resTable))
    if(debug) print(class(resTable[,1]))
    #if(debug) print(resTable)
    if(debug) print(class(resTable[,3]))
    #resTable[resTable[,3]==minValue, 3] = NA
  }
  if (mode=="number") {
    corMat = cor(t(exprMat))
    resTable = melt(corMat)[,1:3]    
  }
  if (mode=="cor") {
    corMat = cor(t(exprMat))
    resTable = melt(corMat)[,1:3]
  }
  if (mode=="p-value") {
    require(Hmisc)
    corRes = rcorr(t(exprMat))
    pvalMat = corRes["P"]
    resTable = melt(pvalMat)[,1:3]
  }
  colnames(resTable) = c("gene1", "gene2", "value")
  
  if (debug) print("sorting table")
  dataTableMode=T
  dataFrameMode=F
  
  if (dataTableMode) {
    resTableNew = as.data.table(resTable)
    resTableNew = resTableNew[order(value, decreasing = T)]
    resTable = as.data.frame(resTableNew)
  }
  if (dataFrameMode) {
    resTable = resTable[order(resTable[,3], decreasing = T),]
  }
  if (!self) {
    resTable = resTable[resTable[,1]!=resTable[,2],]
  }  
  if (debug) print("before cut-off")
  if (debug) print(dim(resTable))
  
  if (debug) print("cut-off data")
  if (mode=="quantile") {
    resTable = resTable[resTable[,3]>=cutOff,]  
  }
  if (mode=="number") {
    resTable = resTable[1:cutOff,]    
  }
  if (mode=="cor") {
    resTable = resTable[resTable[,3]>=cutOff,]
  }
  if (mode=="p-value") {
    resTable = resTable[resTable[,3]<=cutOff,]
  } 
  resTable = resTable[!is.na(resTable[,3]),]
  resTable = resTable[resTable[,3]!=minValue,]
  
  if (debug) print("final")
  if (debug) print(dim(resTable))
  if (debug) print("done")  
  return(resTable)
}


# getClustersList function --> finds communities in a graph (densely connected) via random walks
# u can adjust modularity algorithm between: cluster_louvain & cluster_walktrap
  getClusterList <- function(coexpNet, debug=T) {

  fc = cluster_walktrap(coexpNet, weights = coexpNet[,3])
  if (debug) print("... network clustered")
  #cluster = as.numeric(membership(fc))
  cluster = fc$membership
  geneCluster = data.frame(gene=V(coexpNet)$name,
                           cluster=cluster,
                           stringsAsFactors=F)
  if (debug) print("... cluster sorted")
  geneClusterList = split.data.frame(geneCluster, geneCluster$cluster)
  geneClusterList = lapply(geneClusterList, "[[", "gene")
  geneClusterSizes = do.call(c, lapply(geneClusterList, length))

  print("... done")
  print(paste("... modularity:", as.character(modularity(fc))))
  print(paste("... no. clusters:", as.character(length(geneClusterList))))
  print(paste("... no. genes of max cluster:", as.character(sort(geneClusterSizes,T)[1])))
  return(geneClusterList)
 }


#
# getModulePrunedFrom function ---> keep cluster with min 5 genes 
getModulePrunedFrom <- function(moduleList, cutCluster, debug=T) {
  currGeneClusterLen = unlist(lapply(moduleList, length))
  names(currGeneClusterLen) = names(moduleList)
  if (debug) print(table(currGeneClusterLen>=cutCluster))
  currGeneClusterCut = names(currGeneClusterLen[currGeneClusterLen>=cutCluster])
  moduleList = moduleList[currGeneClusterCut]
  return(moduleList)  
}

#
# getClusterScores - TRANSITIVITY OR CLUSTERING COEFICIENT 
# probability that the adjacent vertices of a vertex are connected
getClusterScores <- function(corNet, geneClusterList) {
  ccScores = lapply(geneClusterList, function(currClusterGenes){
    #currGraph = igraph::subgraph(corNet, currClusterGenes)
    #currGraph = subgraph(corNet, currClusterGenes)
    currGraph = induced.subgraph(corNet, currClusterGenes) # use this one if subgraph is giving you an error
    cc = transitivity(currGraph, "globalundirected")
    return(cc)
  })
  ccScores = unlist(ccScores)
  names(ccScores) = names(geneClusterList)
  return(ccScores)
}


#
# Interaction between modules
areModuleInteract <- function(adjMat, degMat, numEdges, genesA, genesB) {
  observedEdges = observedEdgesBtw(adjMat, genesA, genesB)
  expectedEdges = expectedEdgesBtw(degMat, numEdges, genesA, genesB)
  #diff=(observedEdges-expectedEdges)/expectedEdges
  diff=observedEdges/expectedEdges
  return(diff)
}

getDegMat <- function(igraphNet) {
  #we regard network as undirected
  degs = igraph::degree(igraphNet, mode="in", loops=F)
  #if using undirected graph mode argument is ignored and loops removed
  #degs = igraph::degree(igraphNet, loops= T)
  names(degs) = V(igraphNet)$name
  degMat = as.matrix(degs) %*% t(as.matrix(degs))
  return(degMat)
}
expectedEdgesBtw <- function(degMat, numEdges, genesA, genesB) {
  degMatSel = degMat[genesA, genesB]/(2*numEdges)
  count = sum(degMatSel)
  return(count)
}
observedEdgesBtw <- function(adjMat, genesA, genesB) {
  adjMatSel = adjMat[genesA, genesB]
  count = sum(adjMatSel)
  return(count)
}

#
# SHOWS THE CONECTION BETWEEN THE MODULES
# getModuleMatrix --> need this one for the next one
getModuleMatrix <- function(corNet, modulePruned, debug=F) {
  numEdges = ecount(corNet)
  adjMat = as_adj(corNet)
  degMat = getDegMat(corNet)
  
  numModules = length(modulePruned)
  if (debug) print(numModules)
  moduleMat = matrix(0, numModules, numModules)
  for (i in 1:(numModules-1)) {
    if (debug) print(i)
    for (j in (i+1):numModules) {
      moduleMat[i,j] = areModuleInteract(adjMat, degMat, numEdges, modulePruned[[i]], modulePruned[[j]])
    }
  }  
  moduleMat[lower.tri(moduleMat)] = t(moduleMat)[lower.tri(moduleMat)]
  rownames(moduleMat) = names(modulePruned)
  colnames(moduleMat) = names(modulePruned)
  return(moduleMat)
}


#
# NODE AND EDGE TABle
getModuleNet <- function(prunedList, coexpNet, debug=T) {
  
  moduleCut=0.1
  
  moduleMat = getModuleMatrix(coexpNet, prunedList)
  if (debug) print(moduleMat[1:3,1:3])
  
  if (debug) print(length(prunedList))
  moduleNames = paste("module", names(prunedList), sep=".")
  
  if (debug) print(moduleMat[1:3,1:3])
  moduleMelt = melt(moduleMat)
  moduleMelt = moduleMelt[moduleMelt[,3]>moduleCut,]
  moduleMelt[] = lapply(moduleMelt, as.character)
  moduleMelt[,1] = paste("module", moduleMelt[,1], sep=".")
  moduleMelt[,2] = paste("module", moduleMelt[,2], sep=".")
  
  modulesShown = unique.default(c(moduleMelt[,1], moduleMelt[,2]))
  
  moduleSizes = unlist(lapply(prunedList, length))
  if (debug) print("...size check")
  if (debug) print(length(moduleSizes))
  
  
  moduleAttr = data.frame(node=moduleNames, size=moduleSizes, stringsAsFactors = F)
  
  modulesNotShown = moduleNames[!moduleNames %in% modulesShown]
  return(list(nodeTable = moduleAttr, edgeTable=moduleMelt, nodesNotShown = modulesNotShown ))
}


#
# write THE FILES FOR CYTOSCPAE
write.edge.cytoscape <- function(edgeTable, nodesNotShown, outFile) {
  lenTable=dim(edgeTable)[1]
  lenNodesNotShown = length
  print(head(edgeTable))
  sink(outFile, append=F)
  cat("source")
  cat("\t")
  cat("target")
  cat("\n")
  
  for (i in 1:lenTable) {
    # ints = unlist(edgeTable[i,])
    cat(edgeTable[i,1])
    cat("\t")
    cat(edgeTable[i,2])
    cat("\n")
  }
  for (node in nodesNotShown) {
    cat(node)
    cat("\n")
  }
  sink() 
}

# NODE TABLE
write.node.cytoscape <- function(nodeTable, outFile) {
  write.table(nodeTable, outFile, row.names = F, quote = F, sep="\t")
}