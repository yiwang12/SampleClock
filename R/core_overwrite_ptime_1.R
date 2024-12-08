#' Function to construct pseudotime ordering of samples using either TSCAN or CCA methods
#'
#' Parameters:
#' @param x                 Matrix of features (rows) x samples (columns)
#' @param sample_metadata   Data frame containing sample-level metadata
#' @param sample_name      Column name in sample_metadata identifying samples (default: 'sample') 
#' @param severity_name    Column name in sample_metadata for severity category (default: 'severity')
#' @param severity_num_name Column name in sample_metadata for numerical severity (default: 'sev.level')
#' @param method           Method to use - either 'tscan' or 'cca'
#' @param cluster         Optional pre-defined cluster assignments 
#' @param reduce          Whether to reduce dimensions before analysis (default: FALSE)
#' @param clusternum      Number of clusters for TSCAN (default: 2)
#' @param startcluster    Starting cluster for TSCAN path construction (default: 1)
#' @param modelNames      Model type for mclust clustering (default: 'VVV')
#' @param near_nbr_num    Number of nearest neighbors for graph construction (default: 10)
#
#' Returns:
#' List containing:
#' - For TSCAN: sample_info (ordering + metadata), orthogonal distances, projections, clustering object
#' - For CCA: sample_info (ordering + metadata), slope coefficient, CCA results, CC1/CC2 scores
#'
#' @export
get_sample_pseudotime = function(x, sample_metadata = NULL, 
                                 sample_name = 'sample', severity_name = 'severity', severity_num_name = 'sev.level',
                                 method = c('tscan','cca'),
                                 cluster = NULL, reduce = F, clustermethod = 'mclust', clusternum = 2, startcluster = 1, modelNames = 'VVV',
                                 near_nbr_num = 10){
  method = match.arg(method)
  
  ## arrange sample metadata to have the same order as x
  match_idx = match(colnames(x), sample_metadata[,sample_name])
  sample_metadata = sample_metadata[match_idx, ]
  
  ## tscan
  if (method == 'tscan'){
    sample_cluster = myexprmclust(x, 
                                reduce = reduce,
                                cluster = cluster,
                                clustermethod = clustermethod,
                                clusternum = clusternum,
                                near_nbr_num = near_nbr_num, 
                                modelNames = modelNames)
    sample_metadata$cluster = sample_cluster$clusterid
    #print(str(sample_cluster))
   
    
    TSCAN_result = try(myTSCANorder(sample_cluster,
                              startcluster = startcluster,
                              orderonly = F,
                              listbranch = F))
    if ('try-error' %in% class(TSCAN_result)){
      stop("Failure in pseudotime construction due to inproper sample clustering! Please select another clustering methods or tune the parameters of current clustering method.")
    }
    #print(str(TSCAN_result))
    
    sample_order = TSCAN_result$TSCANorder
    orth_dist = TSCAN_result$orth_dist
    orth_proj = TSCAN_result$orth_proj
    
    
    sample_info = sample_metadata %>%
      inner_join(sample_order, by = setNames("sample_name",sample_name)) %>%
      dplyr::rename(pseudotime = Pseudotime) %>%
      arrange(pseudotime)
      
    return(list(sample_info = sample_info, orth_dist = orth_dist, orth_proj = orth_proj, mclustobj = sample_cluster))
  }
  if (method == 'cca'){
    if (reduce) {
      pc = prcomp(t(x), scale = T)$x[,1:min(5, nrow(x), ncol(x))]###!!!!!!!used top 5 PCs instead of top 2 PCs
    }
    else {
      pc = t(x)
    }
    severity_num = sample_metadata[,severity_num_name]
    print(str(pc))
    print(severity_num)
    res = cancor(pc,severity_num)##???
    #print(str(res))
    #proj = (res$xcoef["PC1",1]*(pc[,"PC1"]-res$xcenter["PC1"]))+ (res$xcoef["PC2",1]*(pc[,"PC2"]-res$xcenter["PC2"])) 
    print(res)
    print(str(pc))
    a = sweep(pc, 2, res$xcenter, '-') #????
    print(str(a))
    proj = sweep(pc, 2, res$xcenter, '-') %*% res$xcoef[,1, drop = F] #?? CC1
    
    CC1 = proj
    CC2 = sweep(pc, 2, res$xcenter, '-') %*% res$xcoef[,2, drop = F] #?? CC1
    
    pseudotime = rank(proj)
    if (cor(pseudotime, severity_num)<0)
      pseudotime = length(pseudotime)+1 -pseudotime
    slope = res$xcoef["PC2",1]/res$xcoef["PC1",1] #???
    #print(paste("slope is", slope))
    #print(paste("corr is ", res$cor))
    
    sample_info = sample_metadata %>%
      mutate(pseudotime = pseudotime) %>%
      arrange(pseudotime)
      
    #print(head(sample_info,10))
    
    return(list(sample_info = sample_info, slope = slope, res=res, CC1=CC1, CC2=CC2))
  }
  
  
}


#' Identify differentially expressed genes along pseudotime trajectory
#'
#' @param x Matrix of expression values (genes x samples)
#' @param sample_metadata Data frame with sample metadata including pseudotime values
#' @param cluster Cell/sample cluster identifier (default: 'all')
#' @param sample_name Column name for sample IDs in metadata (default: 'sample')
#' @param pseudotime_name Column with pseudotime values in metadata (default: 'pseudotime')
#' @param gene_mapping Whether to map gene symbols to IDs using org.Hs.eg.db (default: FALSE)
#' @param deg_parallel Whether to use parallel processing (default: FALSE)
#' @param deg_num_cores Number of cores for parallel processing (default: 10)
#'
#' @return Data frame containing:
#'   - gene: gene identifiers
#'   - pval: raw p-values from GAM model fit
#'   - fdr: adjusted p-values (FDR)
#'   - effect_size: magnitude of expression changes
#'   - cell_cluster: cluster identifier
#'   - signif: boolean indicating significance (FDR < 0.05)
#'   If gene_mapping=TRUE, additional columns:
#'   - ENSEMBL: Ensembl gene IDs
#'   - ENTREZID: Entrez gene IDs
#'
#' @details
#' Uses generalized additive models (GAM) to test for gene expression changes
#' along pseudotime. Expression patterns are modeled using cubic splines.
#' Effect size is calculated as max-min fitted values normalized by residual std dev.
#'
#' @export
run_pseudo_diff_gene = function(x, sample_metadata = NULL, cluster = 'all',
                                sample_name = 'sample', pseudotime_name = 'pseudotime', 
                                gene_mapping = F,
                                deg_parallel = F,
                                deg_num_cores = 10
){
  ptime = sample_metadata[,"sev.level"]
  
  # ptime = sample_metadata[,pseudotime_name]
  sample_metadata = sample_metadata[order(ptime),]
  ptime_order=sample_metadata[,"sev.level"]
  
  genes = rownames(x)
  x_reorder = x[,match(sample_metadata[,sample_name],colnames(x)) %>% na.omit()]
  
  #kp_=order(rowVars(x_reorder))[1:100]
  if(deg_parallel){
    res = mclapply(1:nrow(x_reorder), function(i){
    # res = mclapply(1:nrow(x_reorder), function(i){
      model = mgcv::gam(x_reorder[i,]~s(ptime,k=3))
      
      # model = mgcv::gam(x_reorder[i,]~s(c(1:ncol(x_reorder)),k=3))
      
      smooth = model$fitted.values
      eff = (max(smooth)-min(smooth))/sqrt(sum(residuals(model)^2/model$df.residual))
      pval = pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F) 
      list(smooth = smooth, pval = pval, eff = eff)
    }, mc.cores = deg_num_cores)
  } else{
    res = mclapply(1:nrow(x_reorder), function(i){
    # res = lapply(1:nrow(x_reorder), function(i){
      # model = mgcv::gam(x_reorder[i,]~s(c(1:ncol(x_reorder)),k=3))
      model = mgcv::gam(x_reorder[i,]~s(ptime,k=3))
      smooth = model$fitted.values
      eff = (max(smooth)-min(smooth))/sqrt(sum(residuals(model)^2/model$df.residual))
      pval = pchisq(model$null.deviance - model$deviance, model$df.null - model$df.residual,lower.tail = F) 
      list(smooth = smooth, pval = pval, eff = eff)
    })
  }
  
  pval = sapply(res, function(i) i$pval)
  fdr = p.adjust(pval, method = 'fdr')
  xsmooth = lapply(res, function(i) i$smooth) %>% do.call(rbind, .)
  eff = sapply(res, function(i) i$eff)
  

  result = data.frame(gene = rownames(x_reorder), pval = pval, fdr = fdr, 
                      effect_size = eff) %>%
    mutate(cell_cluster = cluster,
           signif = (fdr < 0.05)) %>%
    arrange(fdr, effect_size)
  
  if (gene_mapping){
    gene_symbols = result$gene
    gene_ids = AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns= c("ENSEMBL", "ENTREZID", "SYMBOL"))
    result = result %>%
      left_join(gene_ids, by = c('gene' = 'SYMBOL'))
  }
  
  return(result)
}

#' Extract features (expression or proportion) for each node in the cell type tree
#'
#' Parameters:
#' @param leaves_info    Data frame containing tree structure information for leaf nodes
#' @param features      Type of features to extract - either 'expr' (expression) or 'prop' (proportion)
#' @param raw_count     Raw count matrix with genes as rows and cells as columns
#' @param cell_meta     Cell metadata containing barcode, sample, and celltype information
#' @param filter_pct    Minimum percentage of samples required (default: 0.9)
#' @param HVG          Whether to identify highly variable genes (default: FALSE)
#' @param qnorm        Whether to quantile normalize data (default: FALSE) 
#' @param batch_correction Whether to perform batch correction (default: FALSE)
#' @param batch        Optional batch information for correction
#' @param clr          Whether to perform centered log-ratio transformation (default: FALSE)
#' @param adj_cov      Additional covariates for adjustment (default: NULL)
#' @param parallel     Whether to use parallel processing (default: TRUE)
#' @param num_cores    Number of cores for parallel processing (default: 10)
#'
#' Returns:
#' For expression features:
#' - hvg: List of highly variable genes per node
#' - agg: Aggregated features across tree levels
#' - all: All features for each node
#'
#' For proportion features:
#' - List of cell type proportions for each tree level
#'
#' @export
get_tree_node_feature = function(leaves_info, features = c('expr','prop'), 
                                 raw_count, cell_meta, filter_pct = 0.9, 
                                 HVG = F, qnorm = F, batch_correction = F, batch = NULL, 
                                 clr = F, adj_cov = NULL, parallel = T, num_cores = 10
){

  unq_id <- unique(leaves_info$label)

  print(unq_id)

  unq_y <- setdiff(unique(leaves_info$y),max(leaves_info$y)) #remove root
  unq_y_root = unique(leaves_info$y)
  tot_sample_size = length(unique(cell_meta$sample))
  
  if (features == 'expr'){
    if (parallel){
      
      pb.ls.all <- mclapply(unq_id, function(tid){
        print(tid)
        node_info <- leaves_info %>% dplyr::filte(label==tid)
        print(node_info)
        #note: cell_meta must contain these columns: barcode, sample, celltype
        sub_meta <- cell_meta %>% dplyr::filte(celltype %in% node_info$children)
  
        sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]
        
        
        unq_sample_size = length(unique(sub_meta$sample))
        if (unq_sample_size < filter_pct * tot_sample_size){
          print(paste("insufficient sample size for node", tid))
          return(NULL)
        }
        # pseudobulk construction (please fill in)
        pb <- get_sample_pb(s = sub_count, pt = sub_meta$sample, 
                            HVG = HVG, qnorm = qnorm, combat = batch_correction, batch = batch)
        print("str(pb)")
        print(str(pb))
        
        return(pb)
      }, mc.cores = num_cores)
      
    } else{
      
      print(unq_id)
      
      pb.ls.all <- lapply(unq_id,function(tid){
        print(tid)
        node_info <- leaves_info %>% dplyr::filte(label==tid)
        print(node_info)
        # note: cell_meta must contain these columns: barcode, sample, celltype
        sub_meta <- cell_meta %>% dplyr::filte(celltype %in% node_info$children)
        sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]

        
        unq_sample_size = length(unique(sub_meta$sample))
        if (unq_sample_size < filter_pct * tot_sample_size){
          print(paste("insufficient sample size for node", tid))
          return(NULL)
        }
        print((unq_sample_size))
        
        # pseudobulk construction (please fill in)
        pb <- get_sample_pb(s = data.matrix(sub_count), pt = sub_meta$sample, 
                            HVG = HVG, qnorm = qnorm, 
                            combat = batch_correction, batch = batch)

        return(pb)
      })
      
      

    }

    pb = lapply(pb.ls.all, function(x) x$all)
    names(pb) = unq_id

    which_kp=which(sapply(pb, function(xx){!is.null(xx)}))
    # pb <- pb[-which(sapply(pb, is.null))] ##!!!
    pb=pb[which_kp]
    # print(99999)
    # print(names(pb))
    
    pb.ls = lapply(pb.ls.all, function(x) x$hvg)
    names(pb.ls) = unq_id
    
    which_kp=which(sapply(pb.ls, function(xx){!is.null(xx)}))
    
    # pb.ls <- pb.ls[-which(sapply(pb.ls, is.null))] #!!!
    pb.ls=pb.ls[which_kp]

    # 
    pb.ls.agg = lapply(unq_y_root, function(ty){
      label_names = leaves_info %>% dplyr::filte(y == ty) %>% pull(label) %>% unique()
      print(label_names)
      pb.ls.sub = pb.ls[which(names(pb.ls) %in% label_names)]
      pb.ls.agg = lapply(1:length(pb.ls.sub), function(i){
        pb = pb.ls.sub[[i]]
        rownames(pb) = paste0(rownames(pb),':', names(pb.ls.sub)[i])
        return(pb)
      })
      samp <- table(unlist(sapply(pb.ls.agg,colnames)))
      samp <- names(samp)[samp==length(pb.ls.agg)]
      pb.ls.agg <- do.call(rbind,lapply(pb.ls.agg,function(pb) pb[,samp]))
    })

    pb.all = list(hvg = pb.ls, agg = pb.ls.agg, all = pb)
    return(pb.all)
    
  }
  
  if(features == 'prop'){
    cell_clu = cell_meta$celltype
    cell_sample = cell_meta$sample
    
    ctp.ls <- lapply(unq_y, function(ty){
      #print(ty)
      node_info <- leaves_info %>% dplyr::filte(y==ty)
      #print(node_info)
      sub_cell_clu = cell_clu[cell_clu %in% node_info$children]
      sub_cell_sample = cell_sample[cell_clu %in% node_info$children]
      
      match_idx = match(sub_cell_clu, node_info$children)
      sub_cell_clu= node_info$label[match_idx]
      
      ctp = get_ctp(clu = sub_cell_clu, pt = sub_cell_sample, 
                    clr = clr, batch_correction = batch_correction, 
                    batch = batch, adj_cov = NULL)
      return(ctp)
    })
    names(ctp.ls) = unq_y
    return(ctp.ls)
  }
  
}
