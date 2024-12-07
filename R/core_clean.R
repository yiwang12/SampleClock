#' Main function for running pseudotime analysis on single-cell data using expression or proportion features
#'
#' @param raw_count Matrix of raw gene expression counts (genes x cells)
#' @param cell_meta Data frame of cell metadata including cell type annotations
#' @param sample_meta Data frame of sample-level metadata
#' @param leaves_info Tree structure information for cell type hierarchy
#' @param features Type of features to use: 'expr', 'prop', or 'expr_and_prop'
#' @param batch Optional batch information for correction
#' @param pb_filter_pct Proportion of samples required for pseudobulk analysis (default: 0.9)
#' @param pb_HVG Whether to use highly variable genes (default: TRUE)
#' @param pb_qnorm Whether to quantile normalize data (default: TRUE) 
#' @param pb_batch_correction Whether to perform batch correction on pseudobulk (default: FALSE)
#' @param pb_parallel Whether to parallelize pseudobulk computation (default: FALSE)
#' @param num_cores Number of cores for parallel processing (default: 10)
#' @param ctp_clr Whether to CLR transform proportions (default: FALSE)
#' @param ctp_adj_cov Additional covariates for proportion adjustment
#' @param ctp_batch_correction Whether to correct batch effects in proportions (default: FALSE)
#' @param sample_name Column name for sample IDs (default: 'sample')
#' @param severity_name Column name for severity category (default: 'severity')
#' @param severity_num_name Column name for numerical severity (default: 'sev.level')
#' @param ptime_method Method for pseudotime: 'cca' or 'tscan' (default: 'cca')
#' @param clustermethod Clustering method for TSCAN (default: 'louvain')
#' @param cluster_id Pre-defined cluster assignments
#' @param reduce Whether to reduce dimensions (default: FALSE)
#' @param clusternum Number of clusters for TSCAN (default: 2)
#' @param startcluster Starting cluster for TSCAN path (default: 1)
#' @param near_nbr_num Number of neighbors for graph construction (default: 3)
#' @param modelNames Model type for mclust (default: 'VVV')
#' @param pc_scale Whether to scale PCs (default: TRUE)
#' @param eval_method Method for evaluating pseudotime quality (default: 'corr')
#'
#' @return List containing:
#' - pb.ls: Pseudobulk expression data
#' - ctp.ls: Cell type proportion data 
#' - ptime.ls: Pseudotime ordering results
#' - sample_meta_ptime: Sample metadata with pseudotime values
#'
#' @export

runPseudotime = function(raw_count, 
                         cell_meta, 
                         sample_meta, 
                         leaves_info,
                         features = 'expr_and_prop', 
                         batch = NULL, 
                         pb_filter_pct = 0.9, 
                         pb_HVG = T, 
                         pb_qnorm = T, 
                         pb_batch_correction = F, 
                         pb_parallel = F,
                         num_cores = 10, 
                         ctp_clr = F, 
                         ctp_adj_cov = NULL, 
                         ctp_batch_correction = F, 
                         sample_name = 'sample', 
                         severity_name = 'severity', 
                         severity_num_name = 'sev.level',
                         ptime_method = 'cca', 
                         clustermethod = 'louvain', 
                         cluster_id = NULL, 
                         reduce = F, 
                         clusternum = 2, 
                         startcluster = 1, 
                         near_nbr_num = 3, 
                         modelNames = 'VVV',
                         pc_scale = T,
                         eval_method = 'corr'
                         
){
  
  if (features %in% c('expr','expr_and_prop')){
    pb.ls = get_tree_node_feature(
      leaves_info, features = 'expr',
      raw_count = raw_count, cell_meta = cell_meta, filter_pct = pb_filter_pct,
      HVG = pb_HVG, qnorm = pb_qnorm, batch_correction = pb_batch_correction, batch = batch,
      parallel = pb_parallel, num_cores = num_cores
    )
    
    if (ptime_method == 'cca'){
      expr.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'expr',
        pb.ls = pb.ls, sample_metadata = sample_meta, ptime_method = 'cca',
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        pc_scale = pc_scale
      )
    }
    if (ptime_method == 'tscan'){
      expr.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'expr',
        pb.ls = pb.ls, sample_metadata = sample_meta, 
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        ptime_method = 'tscan', clustermethod = clustermethod, modelNames = modelNames,
        cluster = cluster_id, reduce = F, clusternum = clusternum, startcluster = startcluster, near_nbr_num = near_nbr_num, 
        pc_scale = T
      )
    }
    
  }
  
  if (features %in% c('prop','expr_and_prop')){
    ctp.ls = get_tree_node_feature(
      leaves_info = leaves_info, features = 'prop',
      raw_count = NULL, cell_meta = cell_meta, filter_pct = NULL,
      clr = ctp_clr, adj_cov = ctp_adj_cov, batch_correction = ctp_batch_correction, batch = batch
    )
    print(str(ctp.ls))
    
    if (ptime_method == 'cca'){
      prop.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'prop',
        ctp.ls = ctp.ls, sample_metadata = sample_meta, ptime_method = 'cca',
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        pc_scale = pc_scale
      )
    }
    if (ptime_method == 'tscan'){
      prop.ptime = pseudotime_tree_node(
        leaves_info = leaves_info, features = 'prop',
        ctp.ls = ctp.ls, sample_metadata = sample_meta, 
        sample_name =  sample_name, severity_name = severity_name, severity_num_name = severity_num_name,
        ptime_method = 'tscan', clustermethod = clustermethod, modelNames = modelNames,
        cluster = cluster_id, reduce = F, clusternum = clusternum, startcluster = startcluster, near_nbr_num = near_nbr_num, 
        pc_scale = T
      )
    }
  }
  
  if (features == 'expr'){
    ptime.ls = expr.ptime
  } else if (features == 'prop'){
    ptime.ls = prop.ptime
  } else if (features == 'expr_and_prop'){
    ptime.ls = c(expr.ptime, prop.ptime)
  }
  
  #print(str(ptime.ls))
  auto_ptime = auto_pseudotime_sel(ptime.ls, ptime_method = ptime_method, eval_method = eval_method)
  sample_meta_ptime = auto_ptime$opt
  opt_idx = which(names(ptime.ls) == names(auto_ptime$eval)[1])
  ptime.ls = ptime.ls[[opt_idx]]
  
  
  if (features == 'expr'){
    ans = list(
      pb.ls = pb.ls,
      ptime.ls = ptime.ls, 
      sample_meta_ptime = sample_meta_ptime
    )
  } else if (features == 'prop'){
    ans = list(
      ctp.ls = ctp.ls,
      ptime.ls = ptime.ls, 
      sample_meta_ptime = sample_meta_ptime
    )
  } else if (features == 'expr_and_prop'){
    ans = list(
      pb.ls = pb.ls, 
      ctp.ls = ctp.ls,
      ptime.ls = ptime.ls,
      sample_meta_ptime = sample_meta_ptime
    )
  }
}


#' Function to identify differentially expressed genes (DEGs) along pseudotime trajectories
#' for different cell populations or nodes in a hierarchical structure
#'
#' Parameters:
#' @param sample_meta     Data frame containing sample metadata including severity and pseudotime information
#' @param pb.ls          List containing pseudobulk expression data, with 'all' component containing
#                      node-specific expression matrices
#' @param leaves_info    Data frame containing hierarchical structure information with columns:
#                      label, x, y, id, and leaf
#' @param severity_name  Column name in sample_meta for severity information (default: 'severity')
#' @param sample_name    Column name in sample_meta identifying samples (default: 'sample')
#' @param pseudotime_name Column name in sample_meta for pseudotime values (default: 'pseudotime')
#' @param gene_mapping   Logical indicating whether to map genes to additional identifiers (default: FALSE)
#' @param deg_parallel   Logical indicating whether to perform parallel processing (default: FALSE)
#' @param deg_num_cores  Number of cores to use for parallel processing (default: 10)
#'
#' Returns:
#' List containing two components:
#' - tab: List of DEG analysis results for each node/population, including:
#'        gene names, p-values, FDR, effect sizes, and significance calls
#' - num_DEG: Summary data frame with number of DEGs per node in the hierarchy
#'
#' Details:
#' The function performs differential expression analysis using GAM models to identify
#' genes that change significantly along pseudotime. For each node in the hierarchy:
#' 1. Fits GAM model to expression data using pseudotime
#' 2. Calculates effect sizes and significance
#' 3. Adjusts for multiple testing
#' 4. Optionally maps genes to additional identifiers
#' 5. Summarizes results in hierarchy context
#'
#' @export
runPseudoDEG = function(sample_meta, 
                        pb.ls,
                        leaves_info, 
                        severity_name = 'severity',
                        sample_name = 'sample', 
                        pseudotime_name = 'pseudotime', 
                        gene_mapping = F,
                        deg_parallel = F,
                        deg_num_cores = 10)
{
  pb.ls.node = pb.ls$all
  
  
  deg_result = lapply(1:length(pb.ls.node), function(i){
      print(names(pb.ls.node)[i])
      res = run_pseudo_diff_gene(x = pb.ls.node[[i]], sample_metadata = sample_meta, cluster = names(pb.ls.node)[i],
                                 gene_mapping = gene_mapping, deg_parallel = deg_parallel, deg_num_cores = deg_num_cores)
      print(head(res))
      return(res)
  })
  names(deg_result) = names(pb.ls.node)
  
  num_DEGs = sapply(1:length(deg_result), function(i){
    num_DEGs = length(which(deg_result[[i]]$signif))
  })
  names(num_DEGs) = names(deg_result)
  print(num_DEGs)
  num_DEGs_col_name = paste0(severity_name, '.num_DEGs')
  
  
  deg_summary = leaves_info %>%
    dplyr::select(label, x, y ,id, leaf) %>%
    unique() %>%
    mutate(!!num_DEGs_col_name := 0)
  match_idx = match(names(num_DEGs), deg_summary$label)
  deg_summary[match_idx, num_DEGs_col_name] = num_DEGs
  
  return(list(tab = deg_result, num_DEG = deg_summary))
}

#' Performs comprehensive pseudotime analysis including trajectory inference and differential expression
#'
#' Parameters:
#' @param raw_count         Matrix of raw gene expression counts (genes x cells)
#' @param cell_meta         Data frame with cell-level metadata
#' @param sample_meta       Data frame with sample-level metadata
#' @param leaves_info       Data frame containing hierarchical clustering tree information
#' @param features         Type of features to use ('expr_and_prop' for both expression and proportions)
#' @param batch            Optional batch information for correction
#' @param pb_filter_pct    Minimum percentage of samples required for pseudobulk analysis (default: 0.9)
#' @param pb_HVG          Whether to use highly variable genes for pseudobulk (default: TRUE)
#' @param pb_qnorm        Whether to quantile normalize pseudobulk data (default: TRUE)
#' @param pb_batch_correction Whether to perform batch correction on pseudobulk (default: FALSE)
#' @param pb_parallel     Whether to parallelize pseudobulk computation (default: FALSE)
#' @param num_cores       Number of cores for parallel processing (default: 10)
#' @param ctp_clr         Whether to apply CLR transformation to cell type proportions (default: FALSE)
#' @param ctp_adj_cov     Optional covariates to adjust in cell type proportions
#' @param ctp_batch_correction Whether to correct batch effects in proportions (default: FALSE)
#' @param sample_name     Column name for sample IDs (default: 'sample')
#' @param severity_name   Column name for severity category (default: 'severity')
#' @param severity_num_name Column name for numerical severity (default: 'sev.level')
#' @param ptime_method    Method for pseudotime ('cca' or 'tscan', default: 'cca')
#' @param clustermethod   Clustering method for TSCAN (default: 'louvain')
#' @param cluster_id      Optional pre-defined cluster assignments
#' @param reduce         Whether to reduce dimensions before analysis (default: FALSE)
#' @param clusternum     Number of clusters for TSCAN (default: 2)
#' @param startcluster   Starting cluster for TSCAN (default: 1)
#' @param near_nbr_num   Number of nearest neighbors for graph construction (default: 3)
#' @param modelNames     Model type for mclust clustering (default: 'VVV')
#' @param pc_scale       Whether to scale PCs (default: TRUE)
#' @param eval_method    Method to evaluate pseudotime quality (default: 'corr')
#' @param pseudotime_name Column name for pseudotime values (default: 'pseudotime')
#' @param gene_mapping   Whether to map genes to standard identifiers (default: FALSE)
#' @param deg_parallel   Whether to parallelize differential expression (default: FALSE)
#' @param deg_num_cores  Number of cores for parallel DEG analysis (default: 10)
#' @param save_deg       Whether to save DEG results to files (default: FALSE)
#' @param sample_embed_plot Whether to generate sample embedding plots (default: FALSE)
#' @param result_path    Path to save results (default: './result/')
#'
#' Returns:
#' List containing:
#' - sample_meta_ptime: Sample metadata with pseudotime
#' - ptime.ls: Pseudotime analysis results 
#' - pb.ls: Pseudobulk analysis results
#' - deg_table: Differential expression results
#' - deg_summary: Summary of differential expression analysis
#'
#' @export
runPseudoAnalysis = function(raw_count, 
                             cell_meta, 
                             sample_meta, 
                             leaves_info,
                             features = 'expr_and_prop', 
                             batch = NULL, 
                             pb_filter_pct = 0.9, 
                             pb_HVG = T, 
                             pb_qnorm = T, 
                             pb_batch_correction = F, 
                             pb_parallel = F,
                             num_cores = 10, 
                             ctp_clr = F, 
                             ctp_adj_cov = NULL, 
                             ctp_batch_correction = F, 
                             sample_name = 'sample', 
                             severity_name = 'severity', 
                             severity_num_name = 'sev.level',
                             ptime_method = 'cca', 
                             clustermethod = 'louvain', 
                             cluster_id = NULL, 
                             reduce = F, 
                             clusternum = 2, 
                             startcluster = 1, 
                             near_nbr_num = 3, 
                             modelNames = 'VVV',
                             pc_scale = T,
                             eval_method = 'corr',
                             pseudotime_name = 'pseudotime', 
                             gene_mapping = F,
                             deg_parallel = F,
                             deg_num_cores = 10,
                             save_deg = F, 
                             sample_embed_plot = F, 
                             result_path = './result/'
                             
){
  
  ptime.args = as.list(match.call())[-1]
  print(ptime.args)
  last_ptime_args = 'eval_method'
  last_ptime_idx = which(names(ptime.args) == last_ptime_args)
  ptime.args = ptime.args[1:last_ptime_idx]
  
  ptime_result = do.call(runPseudotime, ptime.args)
  
  sample_meta_ptime = ptime_result$sample_meta_ptime$sample_info
  ptime.ls = ptime_result$ptime.ls
  pb.ls = ptime_result$pb.ls
  
  if(sample_embed_plot){
    if (ptime_method == 'cca'){
      p = create_pseudotime_plot(sample_info = sample_meta_ptime, 
                                 rd = ptime.ls$rd, 
                                 title = 'CCA', 
                                 slope = ptime.ls$ptime$slope, 
                                 method = 'cca')
    } else if (ptime_method == 'tscan'){
      p = create_pseudotime_plot(sample_info = sample_meta_ptime, 
                                 rd = ptime.ls$rd, 
                                 mclustobj = ptime.ls$ptime$mclustobj, 
                                 title = 'TSCAN', 
                                 method = 'tscan')
    }
    
    dir.create(result_path)
    dir.create(file.path(result_path,'pseudotime'))
    pdf(file.path(result_path, 'pseudotime', paste0('sample_embed_', ptime_method , '.pdf')), 
        height = 16, width = 16 * ifelse(ptime_method == 'cca', 2, 3))
    print(p)
    dev.off()
    write.csv(sample_meta_ptime,file.path(result_path,'pseudotime', 'pseuodtime.csv'))
  }
  

  deg_result = runPseudoDEG(sample_meta = sample_meta_ptime, pb.ls = pb.ls,
                            leaves_info = leaves_info, severity_name = severity_name,
                            sample_name = sample_name, 
                            pseudotime_name = pseudotime_name, 
                            deg_parallel = deg_parallel,
                            deg_num_cores = deg_num_cores)
  deg_table = deg_result$tab
  deg_summary = deg_result$num_DEG
  
  ans = c(ptime_result, list(deg_table = deg_table,deg_summary = deg_summary))
  
  if(save_deg){
  
    dir.create(result_path)
    dir.create(file.path(result_path, 'deg'))
    
    ## save deg
    for(i in 1:length(deg_table)){
      clu = names(deg_table)[i]
      write.csv(deg_table[[i]], file.path(result_path, 'deg', paste0(clu, '.csv')))
    }
  }
  
  return(ans)
}

#' Function to save and process pseudotime analysis results, including differential expression
#' and sample embeddings visualization
#'
#' Parameters:
#' @param ptime_result      List containing pseudotime analysis results
#' @param leaves_info       Information about the hierarchical structure of cell types
#' @param ptime_method      Method used for pseudotime ('cca' or 'tscan', default: 'cca')
#' @param pseudotime_name   Column name for pseudotime values (default: 'pseudotime')
#' @param sample_name       Column name identifying samples (default: 'sample')
#' @param severity_name     Column name for severity information (default: 'severity')
#' @param gene_mapping      Whether to perform gene ID mapping (default: FALSE)
#' @param deg_parallel      Whether to run differential expression in parallel (default: FALSE)
#' @param deg_num_cores    Number of cores for parallel processing (default: 10)
#' @param save_deg         Whether to save differential expression results (default: FALSE)
#' @param sample_embed_plot Whether to generate sample embedding plots (default: FALSE)
#' @param result_path      Directory path to save results (default: './result/')
#'
#' Returns:
#' List containing:
#' - Original pseudotime results
#' - Differential expression tables (deg_table)
#' - Summary of differential expression results (deg_summary)
#'
#' Details:
#' The function performs several tasks:
#' 1. Extracts pseudotime ordering information
#' 2. Optionally creates and saves visualization plots
#' 3. Runs differential expression analysis
#' 4. Saves results to specified directory structure
#'
#' @export
saveResults = function(ptime_result = NULL,
                       leaves_info = NULL,
                       ptime_method = 'cca',
                       pseudotime_name = 'pseudotime', 
                       sample_name = 'sample',
                       severity_name = 'severity',
                       gene_mapping = F,
                       deg_parallel = F,
                       deg_num_cores = 10,
                       save_deg = F, 
                       sample_embed_plot = F, 
                       result_path = './result/'){
  
  sample_meta_ptime = ptime_result$sample_meta_ptime$sample_info
  ptime.ls = ptime_result$ptime.ls
  pb.ls = ptime_result$pb.ls
  
  if(sample_embed_plot){
    if (ptime_method == 'cca'){
      p = create_pseudotime_plot(sample_info = sample_meta_ptime, 
                                 rd = ptime.ls$rd, 
                                 title = 'CCA', 
                                 slope = ptime.ls$ptime$slope, 
                                 method = 'cca')
    } else if (ptime_method == 'tscan'){
      p = create_pseudotime_plot(sample_info = sample_meta_ptime, 
                                 rd = ptime.ls$rd, 
                                 mclustobj = ptime.ls$ptime$mclustobj, 
                                 title = 'TSCAN', 
                                 method = 'tscan')
    }
    
    dir.create(result_path)
    dir.create(file.path(result_path,'pseudotime'))
    pdf(file.path(result_path, 'pseudotime', paste0('sample_embed_', ptime_method , '.pdf')), 
        height = 16, width = 16 * ifelse(ptime_method == 'cca', 2, 3))
    print(p)
    dev.off()
    write.csv(sample_meta_ptime,file.path(result_path,'pseudotime', 'pseuodtime.csv'))
  }
  
  
  deg_result = runPseudoDEG(sample_meta = sample_meta_ptime, pb.ls = pb.ls,
                            leaves_info = leaves_info, severity_name = severity_name,
                            sample_name = sample_name, 
                            pseudotime_name = pseudotime_name, 
                            deg_parallel = deg_parallel,
                            deg_num_cores = deg_num_cores)
  deg_table = deg_result$tab
  deg_summary = deg_result$num_DEG
  
  ans = c(ptime_result, list(deg_table = deg_table,deg_summary = deg_summary))
  
  if(save_deg){
    
    dir.create(result_path)
    dir.create(file.path(result_path, 'deg'))
    
    ## save deg
    for(i in 1:length(deg_table)){
      clu = names(deg_table)[i]
      write.csv(deg_table[[i]], file.path(result_path, 'deg', paste0(clu, '.csv')))
    }
  }
  
  return(ans)
  
}

