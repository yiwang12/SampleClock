


#' Constructs pseudotime ordering for cell types in a hierarchical tree structure
#'
#' Parameters:
#' @param leaves_info       Data frame with tree structure information (label, y, children columns)
#' @param features         Feature type to analyze - either 'expr' (gene expression) or 'prop' (cell proportions)
#' @param raw_count       Raw count matrix for expression data
#' @param cell_meta       Cell-level metadata including barcode, sample, celltype
#' @param pb.ls           List of pseudobulk data for expression analysis
#' @param ctp.ls          List of cell type proportions for proportion analysis  
#' @param sample_metadata Sample-level metadata
#' @param sample_name     Column name for sample IDs (default: 'sample')
#' @param severity_name   Column name for severity category (default: 'severity')
#' @param severity_num_name Column name for numerical severity (default: 'sev.level')
#' @param ptime_method    Method for pseudotime - 'tscan' or 'cca'
#' @param cluster         Optional pre-defined cluster assignments
#' @param reduce          Whether to reduce dimensions (default: FALSE)
#' @param clusternum      Number of clusters for TSCAN (default: 2)
#' @param clustermethod   Clustering method - 'mclust', 'kmeans', or 'louvain'
#' @param startcluster    Starting cluster for TSCAN (default: 1)
#' @param modelNames      Model type for mclust (default: 'VVV')
#' @param near_nbr_num    Number of nearest neighbors (default: 10)
#' @param pc_scale       Whether to scale PCA (default: FALSE)
#'
#' Returns:
#' List containing pseudotime analysis results for each tree level:
#' - For expression data: PCA coordinates, pseudotime ordering
#' - For proportion data: PCA coordinates, pseudotime ordering 
#'
#' @export

pseudotime_tree_node = function(leaves_info, features = c('expr','prop'), 
                                 raw_count, cell_meta, pb.ls, ctp.ls, 
                                 sample_metadata = NULL,
                                 sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                 ptime_method = c('tscan','cca'), cluster = NULL, reduce = F, clusternum = 2, clustermethod = 'mclust',
                                 startcluster = 1, modelNames = 'VVV', near_nbr_num = 10, pc_scale = F
){
  unq_id <- unique(leaves_info$label)
  unq_y <- setdiff(unique(leaves_info$y),max(leaves_info$y)) #remove root
  unq_y_root = unique(leaves_info$y)
  
  
  if (features == 'expr'){ #' for each cell type, use that cell type's pseudo-bulk data to calculate pseudotime
    pb.ls.agg = pb.ls$agg
    ptime.ls = lapply(1:length(pb.ls.agg), function(i){
      pb = pb.ls.agg[[i]]
      rvar = apply(pb, 1, var)
      pb = pb[which(rvar > 0),] #' 
      pbpca = try(prcomp(t(pb),scale. = pc_scale)$x[,1:min(5, nrow(pb), ncol(pb))])
      if ('try-error' %in% class(pbpca))
        return(NULL)
      pbpca = pbpca[rownames(pbpca) %in% sample_metadata[,sample_name],]
      
      sample_ptime_result = get_sample_pseudotime(x = t(pbpca), sample_metadata = sample_metadata, 
                                                sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                method = ptime_method, cluster = cluster, clustermethod = clustermethod,
                                                reduce = reduce, clusternum = clusternum, startcluster = startcluster, modelNames = modelNames,
                                                near_nbr_num = near_nbr_num)
      return(list(ptime = sample_ptime_result, rd = pbpca))
    })
    names(ptime.ls) = paste0('expr_', unq_y_root)
    return(ptime.ls)
  }
  
  if(features == 'prop'){
    ptime.ls = lapply(1:length(ctp.ls), function(i){

      ctp = ctp.ls[[i]]
      rvar = apply(ctp, 1, var)
      ctp = ctp[which(rvar > 0),]
      ctpca = prcomp(t(ctp),scale. = pc_scale)$x[,1:min(5, nrow(ctp), ncol(ctp))]                 
      ctpca = ctpca[rownames(ctpca) %in% sample_metadata[,sample_name],]
      print(str(ctpca))
      
      sample_ptime_result = get_sample_pseudotime(x = t(ctpca), sample_metadata = sample_metadata, 
                                                sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                method = ptime_method, cluster = cluster, clustermethod = clustermethod,
                                                reduce = reduce, clusternum = clusternum, startcluster = startcluster, modelNames = modelNames,
                                                near_nbr_num = near_nbr_num)
      return(list(ptime = sample_ptime_result, rd = ctpca))
    })
    names(ptime.ls) = paste0('prop_', unq_y)
    return(ptime.ls)
  }
}

#' Calculate pseudotime ordering for hierarchical cell populations using top 2 principal components
#'
#' Parameters:
#' @param leaves_info     Data frame containing tree structure information
#' @param features       Type of features to use - either 'expr' (expression) or 'prop' (proportions)
#' @param raw_count     Raw count matrix for expression data
#' @param cell_meta     Cell-level metadata
#' @param pb.ls         List containing pseudobulk data
#' @param ctp.ls        List containing cell type proportion data
#' @param sample_metadata Data frame with sample-level metadata
#' @param sample_name   Column name for sample IDs (default: 'sample')
#' @param severity_name Column name for severity category (default: 'severity')
#' @param severity_num_name Column name for numerical severity (default: 'sev.level')
#' @param ptime_method  Method for pseudotime - 'tscan' or 'cca'
#' @param cluster      Optional pre-defined cluster assignments
#' @param reduce       Whether to reduce dimensions (default: FALSE)
#' @param clusternum   Number of clusters for TSCAN (default: 2)
#' @param clustermethod Clustering method to use (default: 'mclust')
#' @param startcluster Starting cluster for path (default: 1)
#' @param modelNames   Model type for mclust (default: 'VVV')  
#' @param near_nbr_num Number of nearest neighbors (default: 10)
#' @param pc_scale     Whether to scale PCs (default: FALSE)
#'
#' Returns:
#' List containing pseudotime results for each hierarchical level:
#' - For expression data: pseudotime ordering and PCA coordinates
#' - For proportion data: pseudotime ordering and PCA coordinates
#'
#' @export
pseudotime_tree_node_2pc = function(leaves_info, features = c('expr','prop'), 
                                raw_count, cell_meta, pb.ls, ctp.ls, 
                                sample_metadata = NULL,
                                sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                ptime_method = c('tscan','cca'), cluster = NULL, reduce = F, clusternum = 2, clustermethod = 'mclust',
                                startcluster = 1, modelNames = 'VVV', near_nbr_num = 10, pc_scale = F
){
  unq_id <- unique(leaves_info$label)
  unq_y <- setdiff(unique(leaves_info$y),max(leaves_info$y)) #remove root
  unq_y_root = unique(leaves_info$y)
  
  
  if (features == 'expr'){ #' for each cell type, use that cell type's pseudo-bulk data to calculate pseudotime
    pb.ls.agg = pb.ls$agg
    ptime.ls = lapply(1:length(pb.ls.agg), function(i){
      pb = pb.ls.agg[[i]]
      rvar = apply(pb, 1, var)
      pb = pb[which(rvar > 0),] #' 
      pbpca = try(prcomp(t(pb),scale. = pc_scale)$x[,1:min(2, nrow(pb), ncol(pb))])
      if ('try-error' %in% class(pbpca))
        return(NULL)
      pbpca = pbpca[rownames(pbpca) %in% sample_metadata[,sample_name],]
      
      sample_ptime_result = get_sample_pseudotime(x = t(pbpca), sample_metadata = sample_metadata, 
                                                  sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                  method = ptime_method, cluster = cluster, clustermethod = clustermethod,
                                                  reduce = reduce, clusternum = clusternum, startcluster = startcluster, modelNames = modelNames,
                                                  near_nbr_num = near_nbr_num)
      return(list(ptime = sample_ptime_result, rd = pbpca))
    })
    names(ptime.ls) = paste0('expr_', unq_y_root)
    return(ptime.ls)
  }
  
  if(features == 'prop'){
    ptime.ls = lapply(1:length(ctp.ls), function(i){
      
      ctp = ctp.ls[[i]]
      rvar = apply(ctp, 1, var)
      ctp = ctp[which(rvar > 0),]
      ctpca = prcomp(t(ctp),scale. = pc_scale)$x[,1:min(2, nrow(ctp), ncol(ctp))]                 
      ctpca = ctpca[rownames(ctpca) %in% sample_metadata[,sample_name],]
      print(str(ctpca))
      
      sample_ptime_result = get_sample_pseudotime(x = t(ctpca), sample_metadata = sample_metadata, 
                                                  sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                                  method = ptime_method, cluster = cluster, clustermethod = clustermethod,
                                                  reduce = reduce, clusternum = clusternum, startcluster = startcluster, modelNames = modelNames,
                                                  near_nbr_num = near_nbr_num)
      return(list(ptime = sample_ptime_result, rd = ctpca))
    })
    names(ptime.ls) = paste0('prop_', unq_y)
    return(ptime.ls)
  }
}


#' Automatically select optimal pseudotime ordering based on evaluation metrics
#'
#' Parameters:
#' @param ptime.ls          List of pseudotime results from get_sample_pseudotime
#' @param ptime_method     Method used to generate pseudotime ('tscan' or 'cca')
#' @param sample_name      Column name identifying samples in metadata (default: 'sample')
#' @param pseudotime_name  Column name for pseudotime values (default: 'pseudotime')
#' @param severity_num_name Column name for numerical severity values (default: 'sev.level')
#' @param eval_method      Method to evaluate pseudotime quality:
#'                        - 'corr': Correlation with severity
#'                        - 'pairwise_order': Concordance of pairwise orderings
#'                        - 'R2': R-squared with severity
#'
#' Returns:
#' List containing:
#' - opt: Optimal pseudotime ordering result
#' - eval: Named vector of evaluation metrics for all orderings
#'
#' @export
auto_pseudotime_sel = function(ptime.ls, ptime_method = c('tscan','cca'),
                               sample_name = 'sample', pseudotime_name = 'pseudotime', severity_num_name = 'sev.level',
                               eval_method = c('corr', 'pairwise_order','R2')
                               ){
  eval_method = match.arg(eval_method)
  
  if (ptime_method == 'tscan'){
    result = sapply(1:length(ptime.ls), function(i){
    rd = ptime.ls[[i]]$rd 
    ptime_result = ptime.ls[[i]]$ptime
    ptime = ptime_result$sample_info
    
    orth_dist = ptime_result$orth_dist
    ans = get_eval_metric(sample_metadata = ptime, rd = rd, orth_dist = orth_dist, 
                          sample_name = sample_name, pseudotime_name = pseudotime_name, severity_num_name = severity_num_name,
                          eval_method = eval_method)
  })
  }
  
  if (ptime_method == 'cca'){
    result = sapply(1:length(ptime.ls), function(i){
    ptime = ptime.ls[[i]]$ptime$sample_info
    ans = get_eval_metric(sample_metadata = ptime, 
                          sample_name = sample_name, pseudotime_name = pseudotime_name, severity_num_name = severity_num_name,
                          eval_method = eval_method)
  })
  }
  
  
  optime = ptime.ls[[which.max(result)]]$ptime ###!!!!!
  names(result) = names(ptime.ls)
  result = sort(result, decreasing = T)
  return(list(opt = optime, eval = result))
}

#' Construct pseudotime ordering of cells using TSCAN algorithm
#'
#' @param mclustobj Object containing clustering results from myexprmclust
#' @param startcluster Starting cluster for path construction (optional)
#' @param MSTorder Pre-defined path through MST (optional)
#' @param orderonly Whether to return only cell ordering without additional info
#' @param flip Whether to reverse the path direction
#' @param listbranch Whether to identify and return branch paths
#' @param divide Whether to allow path splitting at non-adjacent clusters
#'
#' @return List containing:
#' \itemize{
#'   \item TSCANorder - Data frame with cell ordering and pseudotime
#'   \item orth_proj - Matrix of cell projections onto path
#'   \item orth_dist - Vector of cell distances from path
#' }
#'
#' @details
#' The function constructs a pseudotime path through clusters by:
#' 1. Finding longest path through MST if not provided
#' 2. Ordering cells between adjacent clusters using projections
#' 3. Optionally identifying branch paths
#' 4. Computing orthogonal distances and projections
#'
#' @export
myTSCANorder <- function(mclustobj,startcluster=NULL,MSTorder = NULL,orderonly=F,flip=F,listbranch=F,divide=T) {
  if (!is.null(MSTorder) & length(MSTorder) == 1) {
    stop("MSTorder is not a path!")
  }
  set.seed(12345)
  clucenter <- mclustobj$clucenter
  row.names(clucenter) <- paste0("clu",1:nrow(clucenter))
  clusterid <- mclustobj$clusterid             
  pcareduceres <- mclustobj$pcareduceres            
  adjmat <- as_adjacency_matrix(mclustobj$MSTtree,sparse=FALSE)
  if (is.null(MSTorder)) {
    orderinMST <- 1
    clutable <- table(mclustobj$clusterid)
    alldeg <- degree(mclustobj$MSTtree)
    if (is.null(startcluster)) {
      allcomb <- expand.grid(as.numeric(names(alldeg)[alldeg==1]),as.numeric(names(alldeg)[alldeg==1]))
      allcomb <- allcomb[allcomb[,1] < allcomb[,2],]
      numres <- t(apply(allcomb, 1, function(i) {
        tmp <- as.vector(get.shortest.paths(mclustobj$MSTtree,i[1],i[2])$vpath[[1]])
        c(length(tmp), sum(clutable[tmp]))
      }))
      optcomb <- allcomb[order(numres[,1],numres[,2],decreasing = T)[1],]
      branchcomb <- allcomb[-order(numres[,1],numres[,2],decreasing = T)[1],,drop=F]
      MSTorder <- get.shortest.paths(mclustobj$MSTtree,optcomb[1],optcomb[2])$vpath[[1]]       
    } else {
      allcomb <- cbind(startcluster,setdiff(as.numeric(names(alldeg)[alldeg==1]),startcluster))
      numres <- t(apply(allcomb, 1, function(i) {
        tmp <- as.vector(get.shortest.paths(mclustobj$MSTtree,i[1],i[2])$vpath[[1]])
        c(length(tmp), sum(clutable[tmp]))
      }))
      optcomb <- allcomb[order(numres[,1],numres[,2],decreasing = T)[1],]
      branchcomb <- allcomb[-order(numres[,1],numres[,2],decreasing = T)[1],,drop=F]
      MSTorder <- get.shortest.paths(mclustobj$MSTtree,optcomb[1],optcomb[2])$vpath[[1]] 
    }
    if (flip) MSTorder <- rev(MSTorder)
  } else {
    edgeinMST <- sapply(1:(length(MSTorder)-1),function(i) {
      adjmat[MSTorder[i],MSTorder[i+1]]
    })
    if (divide) {
      if (sum(edgeinMST==0) > 0) {
        orderinMST <- 0
      } else {
        orderinMST <- 1
      }      
    } else {
      orderinMST <- 0
    }
  }          
  
  
  internalorderfunc <- function(internalorder,MSTinout) {
    TSCANorder <- NULL
    orth_proj<- NULL
    orth_dist<- NULL
    
    for (i in 1:(length(internalorder)-1)) {                  
      currentcluid <- internalorder[i]
      nextcluid <- internalorder[i + 1]
      currentclucenter <- clucenter[currentcluid,]
      nextclucenter <- clucenter[nextcluid,]
      
      currentreduceres <- pcareduceres[clusterid==currentcluid,]
      if (MSTinout) {
        connectcluid <- as.numeric(names(which(adjmat[currentcluid,] == 1)))      
      } else {
        if (i == 1) {
          connectcluid <- nextcluid      
        } else {
          connectcluid <- c(nextcluid,internalorder[i - 1])
        }                        
      }
      
      cludist <- sapply(connectcluid, function(x) {                              
        rowSums(sweep(currentreduceres,2,clucenter[x,],"-")^2)
      })
      mindistid <- apply(cludist,1,which.min)
      
      
      edgecell <- names(which(mindistid == which(connectcluid == nextcluid)))
     
      
      
      difvec <- nextclucenter - currentclucenter
      tmppos <- pcareduceres[edgecell,,drop=F] %*% difvec
      pos <- as.vector(tmppos)
      names(pos) <- row.names(tmppos)
     
      
      P <- difvec %*% t(difvec)/as.numeric(t(difvec) %*% difvec)
      edgecell_coord <- pcareduceres[edgecell,,drop=F]
      print(str(edgecell_coord))
      orth_proj_tmp <- t(sapply(1:nrow(edgecell_coord), function(r){
        P %*% edgecell_coord[r, ] + (diag(nrow(P))-P) %*% currentclucenter
      }))
      rownames(orth_proj_tmp) = rownames(edgecell_coord) 
      orth_dist_tmp = sapply(1:nrow(orth_proj_tmp), function(r){
        sum((edgecell_coord[r,]-orth_proj_tmp[r,])^2)
      })
      names(orth_dist_tmp) = rownames(edgecell_coord)

      orth_proj = rbind(orth_proj, orth_proj_tmp)
      orth_dist = c(orth_dist, orth_dist_tmp)
          
      
      
      TSCANorder <- c(TSCANorder,names(sort(pos)))  
      
      nextreduceres <- pcareduceres[clusterid==nextcluid,,drop=F]     
      if (MSTinout) {
        connectcluid <- as.numeric(names(which(adjmat[nextcluid,] == 1)))
      } else {
        if (i == length(internalorder)-1) {
          connectcluid <- currentcluid      
        } else {
          connectcluid <- c(currentcluid,internalorder[i + 2])
        }                        
      }
      
      cludist <- sapply(connectcluid, function(x) { 
        rowSums(sweep(nextreduceres,2,clucenter[x,],"-")^2)
      })
      if (length(cludist)==1) {
        mindistid <- 1
      } else {
        mindistid <- apply(cludist,1,which.min)
      }
      
      edgecell <- names(which(mindistid == which(connectcluid == currentcluid)))
      
      difvec <- nextclucenter - currentclucenter
      tmppos <- pcareduceres[edgecell,,drop=F] %*% difvec
      pos <- as.vector(tmppos)
      names(pos) <- row.names(tmppos)
      
      
      P <- difvec %*% t(difvec)/as.numeric(t(difvec) %*% difvec)
      edgecell_coord <- pcareduceres[edgecell,,drop=F]
      orth_proj_tmp <- t(sapply(1:nrow(edgecell_coord), function(r){
        P %*% edgecell_coord[r, ] + (diag(nrow(P))-P) %*% currentclucenter
      }))
      rownames(orth_proj_tmp) = rownames(edgecell_coord) 
      
      orth_dist_tmp = sapply(1:nrow(orth_proj_tmp), function(r){
        sum((edgecell_coord[r,]-orth_proj_tmp[r,])^2)
      })
      names(orth_dist_tmp) = rownames(edgecell_coord)
      
      orth_proj = rbind(orth_proj, orth_proj_tmp)
      orth_dist = c(orth_dist, orth_dist_tmp)
      
      TSCANorder <- c(TSCANorder,names(sort(pos)))
     
      
    }
    
    if (orderonly) {
      TSCANorder = TSCANorder      
    } else {
      TSCANorder = data.frame(sample_name=TSCANorder,State=clusterid[TSCANorder],Pseudotime=1:length(TSCANorder),stringsAsFactors = F)            
    }
    
    return(list(TSCANorder = TSCANorder, orth_proj = orth_proj, orth_dist = orth_dist))
  }
  if (!orderinMST) {
    internalorderfunc(MSTorder,0)            
  } else {
    if (exists("branchcomb") & listbranch) {                  
      allres <- list()
      backbone_result = internalorderfunc(MSTorder,1)
      allres["TSCANorder"] = list();allres[["orth_proj"]] = list(); allres[["orth_dist"]] = list()
      
      allres[["TSCAN_order"]][[paste("backbone",paste(MSTorder,collapse = ','))]] <- backbone_result$TSCANorder
      allres[["orth_proj"]][[paste("backbone",paste(MSTorder,collapse = ','))]] = backbone_result$orth_proj
      allres[["orth_dist"]][[paste("backbone",paste(MSTorder,collapse = ','))]] = backbone_result$orth_dist
      
      for (tmpcombid in 1:nrow(branchcomb)) {
        tmporder <- get.shortest.paths(mclustobj$MSTtree,branchcomb[tmpcombid,1],branchcomb[tmpcombid,2])$vpath[[1]] 
        if (flip) tmporder <- rev(tmporder)
        branch_result = internalorderfunc(tmporder,1)
        allres[["TSCAN_order"]][[paste("branch:",paste(tmporder,collapse = ','))]] <- branch_result$TSCANorder 
        allres[["orth_proj"]][[paste("branch:",paste(tmporder,collapse = ','))]] <- branch_result$orth_proj
        allres[["orth_dist"]][[paste("branch:",paste(tmporder,collapse = ','))]] <- branch_result$orth_dist
      }
      allres
    } else {
      internalorderfunc(MSTorder,1)            
    }      
  }
}

