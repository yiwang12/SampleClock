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



# `x`: the data matrix of features $\times$ samples, i.e. each row is a feature of gene expression or cell type proportion, each column is a patient sample.

# `sample_metadata`: patient-level metadata, including several important variables.

# `sample_name` and `severity_name`: the variable name in `sample_metadata` that indicates the sample and severity. The default setting of `sample_name` is "Patient" and `severity_name` is "severity". 

# `method`: users now have two options, either using TSCAN for sample-level pseudotime construction or use CCA for optimal severity axis detection. 

# `cluster = NULL, reduce = F, clusternum = 2, startcluster = 1`: parameters for TSCAN.

# `severity_num_name`: the variable name in 'sample_metadata' that indicates the numerical valaue of severity severity. 

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

# pb.ls = get_tree_node_feature(
#   leaves_info
#   features = 'expr'
#   raw_count = raw_count
#   cell_meta = cell_meta
#   filter_pct = pb_filter_pc
#   HVG = pb_HVG
#   qnorm = pb_qnorm
#   batch_correction = pb_batch_correction
#   batch = batch
#   parallel = pb_parallel
#   num_cores = num_cores
# )
# filter_pct = 0.9 
# HVG = F
# qnorm = F
# batch_correction = F
# batch = NULL
# clr = F
# adj_cov = NULL
# parallel = T
# num_cores = 10
# raw_count=raw_count_sub
# cell_meta=cell_meta_sub
get_tree_node_feature = function(leaves_info, features = c('expr','prop'), 
                                 raw_count, cell_meta, filter_pct = 0.9, 
                                 HVG = F, qnorm = F, batch_correction = F, batch = NULL, 
                                 clr = F, adj_cov = NULL, parallel = T, num_cores = 10
){
  print(111111)
  
  unq_id <- unique(leaves_info$label)
  print(2222)
  
  unq_y <- setdiff(unique(leaves_info$y),max(leaves_info$y)) #remove root
  print(3333)
  unq_y_root = unique(leaves_info$y)
  print(4444)
  tot_sample_size = length(unique(cell_meta$sample))
  
  
  if (features == 'expr'){
    if (parallel){
      
      pb.ls.all <- mclapply(unq_id, function(tid){
        print(tid)
        node_info <- leaves_info %>% filter(label==tid)
        print(node_info)
        # note: cell_meta must contain these columns: barcode, sample, celltype
        sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
  
        sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]
        
        
        unq_sample_size = length(unique(sub_meta$sample))
        if (unq_sample_size < filter_pct * tot_sample_size){
          print(paste("insufficient sample size for node", tid))
          return(NULL)
        }
        # pseudobulk construction (please fill in)
        pb <- get_sample_pb(s = sub_count, pt = sub_meta$sample, 
                            HVG = HVG, qnorm = qnorm, combat = batch_correction, batch = batch)
        print(str(pb))
        return(pb)
      }, mc.cores = num_cores)
      
    } else{
      
      print(5555)
      print(unq_id)
      
      pb.ls.all <- lapply(unq_id,function(tid){
        print(tid)
        node_info <- leaves_info %>% filter(label==tid)
        print(node_info)
        # note: cell_meta must contain these columns: barcode, sample, celltype
        sub_meta <- cell_meta %>% filter(celltype %in% node_info$children)
        sub_count <- raw_count[, match(sub_meta$barcode, colnames(raw_count))]
        # print("head(sub_meta)")
        # print(head(sub_meta))
        print("00000")
        print(dim(sub_meta))
        print(dim(sub_count))
        print(head(sub_meta))
        print(head(sub_count[1:2,1:2]))
        
        print("00000")
        # print("head(match(sub_meta$barcode, colnames(raw_count)))")
        # print(head(match(sub_meta$barcode, colnames(raw_count))))
        
        # print((unq_sample_size))
        
        unq_sample_size = length(unique(sub_meta$sample))
        if (unq_sample_size < filter_pct * tot_sample_size){
          print(paste("insufficient sample size for node", tid))
          return(NULL)
        }
        print((unq_sample_size))
        
        # pseudobulk construction (please fill in)
        pb <- get_sample_pb(s = sub_count, pt = sub_meta$sample, 
                            HVG = HVG, qnorm = qnorm, 
                            combat = batch_correction, batch = batch)
        print(str(pb))
        # 
        # print("00001110")
        # print(names(pb[[1]]))
        # ls(pb)
        print(head(pb$all))
        return(pb)
      })
      
      
      print(66666)
      
    }
    print(77777)
    
    print(names(pb.ls.all))
    # NULL
    pb = lapply(pb.ls.all, function(x) x$all)
    names(pb) = unq_id
    print(88888)
    
    print(names(pb))
    
    
    which_kp=which(sapply(pb, function(xx){!is.null(xx)}))
    # pb <- pb[-which(sapply(pb, is.null))] ##!!!
    pb=pb[which_kp]
    print(99999)
    print(names(pb))
    
    pb.ls = lapply(pb.ls.all, function(x) x$hvg)
    names(pb.ls) = unq_id
    
    which_kp=which(sapply(pb.ls, function(xx){!is.null(xx)}))
    
    # pb.ls <- pb.ls[-which(sapply(pb.ls, is.null))] #!!!
    pb.ls=pb.ls[which_kp]
    
    print(3333333)
    print(unq_y_root)
    # [1] 2 1 0
    # print(label)
    print(names(pb.ls))
    print(3333333)
    print(3333333)
    
    pb.ls.agg = lapply(unq_y_root, function(ty){
      label_names = leaves_info %>% filter(y == ty) %>% pull(label) %>% unique()
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
    # print(44444444)
    # 
    # names(pb.ls.agg) = unq_y_root
    # print(str(pb.ls.agg))
    # pb.ls.agg=NA
    pb.all = list(hvg = pb.ls, agg = pb.ls.agg, all = pb)
    return(pb.all)
    
  }
  
  if(features == 'prop'){
    cell_clu = cell_meta$celltype
    cell_sample = cell_meta$sample
    
    ctp.ls <- lapply(unq_y, function(ty){
      #print(ty)
      node_info <- leaves_info %>% filter(y==ty)
      #print(node_info)
      sub_cell_clu = cell_clu[cell_clu %in% node_info$children]
      sub_cell_sample = cell_sample[cell_clu %in% node_info$children]
      
      match_idx = match(sub_cell_clu, node_info$children)
      sub_cell_clu= node_info$label[match_idx]
      
      ctp = get_ctp(clu = sub_cell_clu, pt = sub_cell_sample, 
                    clr = clr, batch_correction = batch_correction, 
                    batch = batch, adj_cov = NULL)
      #print(str(ctp))
      return(ctp)
    })
    names(ctp.ls) = unq_y
    return(ctp.ls)
  }
  
}

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
  
  
  if (features == 'expr'){ # for each cell type, use that cell type's pseudo-bulk data to calculate pseudotime
    pb.ls.agg = pb.ls$agg
    ptime.ls = lapply(1:length(pb.ls.agg), function(i){
      pb = pb.ls.agg[[i]]
      rvar = apply(pb, 1, var)
      pb = pb[which(rvar > 0),] # 
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
  
  
  if (features == 'expr'){ # for each cell type, use that cell type's pseudo-bulk data to calculate pseudotime
    pb.ls.agg = pb.ls$agg
    ptime.ls = lapply(1:length(pb.ls.agg), function(i){
      pb = pb.ls.agg[[i]]
      rvar = apply(pb, 1, var)
      pb = pb[which(rvar > 0),] # 
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