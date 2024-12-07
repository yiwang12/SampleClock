library(TSCAN)
library(TreeCorTreat)
library(tidyverse)
library(Matrix)
library(sva)
library(preprocessCore)
library(mgcv)
library(org.Hs.eg.db)
library(parallel)
library(igraph)
library(mclust)
library(scran)

# Performs clustering on expression data using various methods and constructs a minimum spanning tree
#
# Parameters:
# @param data           Expression matrix with features as rows and samples as columns
# @param clustermethod  Clustering method to use: 'mclust', 'kmeans', or 'louvain' (default: 'mclust')
# @param clusternum     Range of cluster numbers to try for mclust, or specific number for kmeans (default: 2:9)
# @param modelNames     Model type for mclust clustering (default: 'VVV')
# @param reduce         Whether to reduce dimensions using PCA before clustering (default: TRUE)
# @param cluster       Optional pre-defined cluster assignments
# @param near_nbr_num  Number of nearest neighbors for Louvain clustering (default: 10)
#
# Returns:
# List containing:
# - pcareduceres: PCA-reduced data (if reduce=TRUE) or original data
# - MSTtree: Minimum spanning tree connecting cluster centers
# - clusterid: Cluster assignments for each sample
# - clucenter: Coordinates of cluster centers
#
# @export
myexprmclust <- function (data, clustermethod = 'mclust', clusternum = 2:9, modelNames = "VVV", reduce = T, cluster = NULL, near_nbr_num = 10) {
  set.seed(12345)
  if (reduce) {
    sdev <- prcomp(t(data), scale = T)$sdev[1:20]
    x <- 1:20
    optpoint <- which.min(sapply(2:10, function(i) {
      x2 <- pmax(0, x - i)
      sum(lm(sdev ~ x + x2)$residuals^2)
    }))
    pcadim = optpoint + 1
    tmpdata <- t(apply(data, 1, scale))
    colnames(tmpdata) <- colnames(data)
    tmppc <- prcomp(t(tmpdata), scale = T)
    pcareduceres <- t(tmpdata) %*% tmppc$rotation[, 1:pcadim]
  }
  else {
    pcareduceres <- t(data)
  }
  print(str(pcareduceres))
  
  
  if (clustermethod=='mclust') {
    if (is.null(cluster)) {   
      clusternum <- clusternum[clusternum > 1]
      res <- suppressWarnings(Mclust(pcareduceres, G = clusternum, modelNames = modelNames))
      clusterid <- apply(res$z, 1, which.max)
      clunum <- res$G
    } else {
      clunum <- length(unique(cluster))
      clusterid <- cluster
    }
  } else if (clustermethod == 'kmeans'){
    clusterid <- kmeans(pcareduceres,clusternum)$cluster
    print(clusterid)
    clunum <- clusternum
  } else if (clustermethod == 'louvain'){
    print(str(pcareduceres))
    print(near_nbr_num)
    snngraph = scran::buildSNNGraph(pcareduceres, transposed= T, k = near_nbr_num, d = NA)
    clusterid = cluster_louvain(snngraph)$membership
    names(clusterid) = rownames(pcareduceres)
    print(clusterid)
    clunum = length(unique(clusterid))
    print(clunum)
  }
  clucenter <- matrix(0, ncol = ncol(pcareduceres), nrow = clunum)
  for (cid in 1:clunum) {
    clucenter[cid, ] <- colMeans(pcareduceres[names(clusterid[clusterid == cid]), , drop = F])
  }
  
  dp <- as.matrix(dist(clucenter))
  gp <- graph.adjacency(dp, mode = "undirected", weighted = TRUE)
  dp_mst <- minimum.spanning.tree(gp)
  list(pcareduceres = pcareduceres, MSTtree = dp_mst, clusterid = clusterid, clucenter = clucenter)
}