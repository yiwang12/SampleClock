library(tidyverse)
library(patchwork)

severity.color = readRDS("/dcl02/hongkai/data/covid/data/current/palette/severity.rds")
study.color = readRDS("/dcl02/hongkai/data/covid/data/current/palette/study.rds")

theme_paper = function(){
  theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(legend.position = "top", legend.key.size = unit(0.3, "in"),legend.text = element_text(size = 30),legend.title=element_text(size = 30),legend.box = "vertical") + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) +
    theme(axis.text.x = element_text(size=30,color="black"),
          axis.text.y = element_text(size=30,color='black'),
          axis.title.x = element_text(size=40,vjust=-1),
          axis.title.y = element_text(size=40,vjust=1),
          plot.margin=unit(c(1,1,1,1),"cm"))+
    theme(plot.title = element_text(size = 40, hjust = 0.5))
}

## sample embed 
create_pseudotime_plot = function(sample_info, rd, mclustobj, title, slope = NULL, 
                                  method = c('cca','tscan')){
  
  method = match.arg(method)
  
  sample_info = rd %>% 
    as.data.frame() %>%
    mutate(sample = rownames(.)) %>%
    inner_join(sample_info, by = 'sample') 
  
  if (method == 'cca'){
    p1 = ggplot(aes(x = PC1, y = PC2, color = severity, label = sample), data= sample_info)+
      geom_point(size = 5)+
      scale_colour_manual(values = severity.color)+
      #geom_text_repel()+
      geom_abline(slope = slope, intercept = 0, linetype = 'dashed', col = 'red', size = 3)+
      guides(colour = guide_legend(nrow=1,byrow=TRUE)) + 
      labs(color = "Severity") +
      ggtitle("Severity")+
      theme_paper()
    
    p2 = ggplot(aes(x = PC1, y = PC2, color = study, label = sample), data= sample_info)+
      geom_point(size = 5)+
      scale_colour_manual(values = study.color)+
      #geom_text_repel()+
      geom_abline(slope = slope, intercept = 0, linetype = 'dashed', col = 'red', size = 3)+
      guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
      labs(color = "Study") +
      ggtitle("Study")+
      theme_paper()
    
    p_combo = (p1 + p2) + 
      plot_annotation(title = title,
                      theme = theme(plot.title = element_text(size = 40, hjust = 0.5))) 
      
  }
  
  if (method == 'tscan'){
    p1 = TSCAN::plotmclust(mclustobj, cell_point_size = 5) + 
      #geom_text_repel(aes(x = PC1, y = PC2, label = sample), data = sample_info) +
      labs(color = "State") +
      ggtitle("State")+
      theme_paper()
    p2 = TSCAN::plotmclust(mclustobj)+
      #geom_text_repel()+
      geom_point(aes(x = PC1, y = PC2, color = severity), size = 5, data= sample_info)+
      scale_colour_manual(values = severity.color, limits = unique(sample_info$severity))+
      #guides(colour = guide_legend( nrow=2,byrow=TRUE)) + 
      labs(color = 'Severity')+
      ggtitle("Severity")+
      theme_paper()
    p3 = TSCAN::plotmclust(mclustobj) + 
      #scale_colour_manual(values = severity.color)+
      #geom_text_repel()+
      geom_point(aes(x = PC1, y = PC2, color = study), size = 5, data= sample_info)+
      guides(colour = guide_legend(nrow=2,byrow=TRUE)) + 
      scale_colour_manual(values = study.color, limits = unique(sample_info$study))+
      labs(color = 'Study')+
      ggtitle("Study")+
      theme_paper()
    p_combo = (p1 + p2 + p3) + 
      plot_annotation(title = title,
                      theme = theme(plot.title = element_text(size = 40, hjust = 0.5))) 
      
  }
  
  return(p_combo)
}


## cell prop
get_cell_prop_traj = function(cell_type_prop, 
                              sample_metadata = NULL,
                              sample_name = 'sample', pseudotime_name = 'pseudotime'){
  
  rn = rownames(cell_type_prop)
  celltype_list = gsub('(.*)_(.*)','\\1', rn)
  cluster_list = gsub('(.*)_(.*)','\\2', rn)
  celltype_label = paste0(cluster_list, '\n', celltype_list)
  patient_order = sample_metadata[, sample_name]
  
  cp = cell_type_prop[, match(patient_order,colnames(cell_type_prop)),drop = F]
  
  cp_df = cp %>%
    as.data.frame() %>%
    mutate(celltype = celltype_label) %>%
    gather(!!sample_name, value, -celltype) 
  str(cp_df)
  
  a = sample_metadata %>% dplyr::select({{sample_name}}, {{pseudotime_name}})
  str(a)
 
  cp_df = cp_df  %>%
    inner_join(sample_metadata %>% dplyr::select({{sample_name}}, {{pseudotime_name}}), by = sample_name)
  
  str(cp_df)
  
  p = cp_df %>%
    mutate(celltype = factor(celltype, levels = celltype_label)) %>%
    ggplot(aes(x = get(pseudotime_name), y = value))+
    geom_point()+ geom_smooth()+ 
    facet_wrap(vars(celltype),scales = "free_y")+
    ylab("Cell Type Proportion")+
    theme_classic()
  
  return(p)
    
}

#topclu = c(1,2,3,4,8,13,14,24,25); names(topclu) = paste0("c",topclu)
#topct = c(rep('T',5), rep('Mono',2), rep('B',2))

#ctp = get_ctp(clu = cell_meta$celltype, pt = cell_meta$sample)
#p1 = get_cell_prop_traj(cell_type_prop = ctp, 
#                        sample_metadata = sample_meta_ptime)
#p1

## gene expr
get_gene_expr_traj = function(x, gene_list, cluster, sample_metadata, 
                              pseudotime_name = 'pseudotime',
                              sample_name = 'sample'){
  ptime = sample_metadata[,pseudotime_name]
  sample_metadata = sample_metadata[order(ptime),]
  match_id = match(sample_metadata[, sample_name], colnames(x))
  x = x[gene_list, match_id , drop = F]
  colnames(x) = 1:ncol(x)
  
 
  df = x %>%
    as.data.frame() %>%
    mutate(gene=rownames(.)) %>%
    pivot_longer(!gene, names_to = 'Pseudotime', values_to = 'Expression') %>%
    mutate(Pseudotime = as.numeric(Pseudotime))

  p = df%>%
    ggplot(aes(x = Pseudotime, y = Expression))+
    geom_point()+ geom_smooth()+
    facet_wrap(vars(gene), scales = "free_y") + 
    ggtitle(cluster)+
    theme_classic()
  return(p)
}
#p = get_gene_expr_traj(x = pb.ls.node[['Mono']], 
#                       gene_list = head(deg_table[['Mono']]$gene), 
#                       cluster = 'Mono', 
#                       sample_metadata = sample_meta_ptime)

