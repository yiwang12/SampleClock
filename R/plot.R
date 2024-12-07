
#' Creates visualization plots for pseudotime analysis using either CCA or TSCAN methods
#'
#' Parameters:
#' @param sample_info    Data frame containing sample metadata and pseudotime information
#' @param rd            Reduced dimension coordinates (typically PCA or other dimensionality reduction)
#' @param mclustobj     Clustering object from TSCAN analysis (required for TSCAN method)
#' @param title         Main title for the combined plot
#' @param slope         Slope for CCA projection line (required for CCA method)
#' @param method        Visualization method to use - either 'cca' or 'tscan'
#'
#' Details:
#' For CCA method:
#' - Creates two plots showing PC1 vs PC2:
#'   1. Points colored by disease severity
#'   2. Points colored by study origin
#' - Includes projection line with specified slope
#'
#' For TSCAN method:
#' - Creates three plots:
#'   1. Points colored by cluster state
#'   2. Points colored by disease severity
#'   3. Points colored by study origin
#' - Includes minimum spanning tree connecting cluster centers
#'
#' All plots use consistent color schemes:
#' - Severity: Green (Healthy) -> Blue (Mild) -> Orange (Moderate) -> Red (Severe)
#' - Studies: Preset color palette for 12 different studies
#'
#' Returns:
#' A combined plot arranged using patchwork with shared title
#'
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import TSCAN
#'
#' @export
create_pseudotime_plot = function(sample_info, rd, mclustobj, title, slope = NULL, 
                                  method = c('cca','tscan')){
  study_names <- c("Aruna", "Guo", "Lee", "Mudd", "Silvin", 
                   "SS1", "SS2", "Su", "Wen", "Wilk", "Yu", "Zhu")
  
  color_codes <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                   "#FFD92F", "#E5C494")
  study.color <- setNames(color_codes, study_names)
  
  
  severity.color = c("#4DAF4A" ,"#377EB8",  "orange", "#E41A1C" )
  names(severity.color)[1:4]=c("Healthy","Mild","Moderate","Severe")
  
  #' Create named vector
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

#' Creates trajectory plots for cell type proportions across pseudotime
#'
#' Parameters:
#' @param cell_type_prop   Matrix of cell type proportions (rows: cell_types, columns: samples)
#' @param sample_metadata  Data frame containing sample metadata
#' @param sample_name     Column name in sample_metadata for sample IDs (default: 'sample')
#' @param pseudotime_name Column name in sample_metadata for pseudotime values (default: 'pseudotime')
#'
#' Returns:
#' ggplot object showing cell type proportion trajectories faceted by cell type
#'
#' @importFrom ggplot2 ggplot geom_point geom_smooth facet_wrap theme_classic
#' @importFrom dplyr select inner_join
#' @importFrom tidyr gather
#' @export
get_cell_prop_traj = function(cell_type_prop, 
                              sample_metadata = NULL,
                              sample_name = 'sample', pseudotime_name = 'pseudotime'){
  study_names <- c("Aruna", "Guo", "Lee", "Mudd", "Silvin", 
                   "SS1", "SS2", "Su", "Wen", "Wilk", "Yu", "Zhu")
  
  color_codes <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                   "#FFD92F", "#E5C494")
  study.color <- setNames(color_codes, study_names)
  
  
  severity.color = c("#4DAF4A" ,"#377EB8",  "orange", "#E41A1C" )
  names(severity.color)[1:4]=c("Healthy","Mild","Moderate","Severe")
  
  # Create named vector
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


#' Creates trajectory plots for gene expression across pseudotime
#'
#' @param x                Matrix of gene expression values (genes x samples)
#' @param gene_list       Vector of gene names to plot
#' @param cluster         Cluster identifier for plot title
#' @param sample_metadata Data frame containing sample metadata
#' @param pseudotime_name Column name in sample_metadata for pseudotime values (default: 'pseudotime')
#' @param sample_name     Column name in sample_metadata identifying samples (default: 'sample')
#'
#' @return ggplot object showing gene expression trajectories across pseudotime
#'   - Each gene gets its own facet
#'   - Points show expression values
#'   - Smoothed trend line added for each gene
#'   - Free y-axis scaling per facet
#'   - Classic theme with customized formatting
#'
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth facet_wrap theme_classic
#' @importFrom dplyr mutate
#'
#' @export
get_gene_expr_traj = function(x, gene_list, cluster, sample_metadata, 
                              pseudotime_name = 'pseudotime',
                              sample_name = 'sample'){
  study_names <- c("Aruna", "Guo", "Lee", "Mudd", "Silvin", 
                   "SS1", "SS2", "Su", "Wen", "Wilk", "Yu", "Zhu")
  
  color_codes <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                   "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
                   "#FFD92F", "#E5C494")
  study.color <- setNames(color_codes, study_names)
  
  
  severity.color = c("#4DAF4A" ,"#377EB8",  "orange", "#E41A1C" )
  names(severity.color)[1:4]=c("Healthy","Mild","Moderate","Severe")
  
  # Create named vector
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
