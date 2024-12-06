## input: 
# raw count 
# cell meta
# sample meta
# tree structure

## output:
# pb.ls, ctp.ls
# sample meta
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

## input: 
# sample meta
# leaves_info: tree structure

## output:
# deg result
# deg summary

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

## unified function
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

