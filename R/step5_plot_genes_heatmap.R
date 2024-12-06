library(pheatmap)
library(matrixStats)
library(ggplot2)
library(cowplot)
library(stringr)

my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
                        "#E16462FF", "#FCA636FF", "#F0F921FF")
                        
###############################################################



#################

plot_heatmap<-function(list_norm, genes_, clu_genes_, pseudotime_){

  pseudotime__ = pseudotime_[colnames(list_norm[[1]])]
  
  samples_kp = which(!is.na(pseudotime__))
  pseudotime__ = pseudotime__[samples_kp]
  
  for(cellType in c(names(list_norm))){
    mat_tmp = list_norm[[cellType]][genes_, samples_kp]
    # mat_tmp=mat_tmp[!is.nan(rowSums(mat_tmp)),]
    # mat_tmp=mat_tmp[,!is.nan(colSums(mat_tmp))]
    
    ####
    mat_tmp_plot=mat_tmp[order(clu_genes_[genes_]),
                         order(pseudotime__, decreasing = F)]
    mat_tmp_plot=t(scale(t(mat_tmp_plot)))
    mat_tmp_plot[mat_tmp_plot>1] = 1
    mat_tmp_plot[mat_tmp_plot< -1] = -1
    
    ### annotate rows by cluster number
    clu_genes_kp = clu_genes_[genes_]
    clu_genes_ordered  = clu_genes_kp[order(clu_genes_kp)]
    annotation_row = data.frame(cluster = as.factor(clu_genes_ordered))
    rownames(annotation_row)=rownames(mat_tmp_plot)
    
    print(colnames(mat_tmp_plot))
    
    ###
    # if(cellType == "M"){cellType_name = "Monocyte"}else{
      # cellType_name = paste0(cellType," cell")
    # }
    pheatmap(mat_tmp_plot,
             cluster_rows=F,cluster_cols=F,
             main = cellType,
             annotation_row = annotation_row)
  }
  
}



plot_consist_div <- function(list_norm, clu_genes, clusters_divergent, clusters_consistent, pseudotime_rna, dir_plot=NULL, interactive=F){
  dir.create(dir_plot)
  
  if(!is.null(dir_plot)){
    pdf(paste0(dir_plot, "/geneCluster_heatmap_div.pdf"))
    genes_tmp = names(clu_genes)[as.character(clu_genes)%in%as.character(clusters_divergent)]
    p1 = plot_heatmap(list_norm, genes_tmp, clu_genes, pseudotime_rna)
    print(p1)
    dev.off()
    
    
    pdf(paste0(dir_plot, "/geneCluster_heatmap_consistent.pdf"))
    genes_tmp = names(clu_genes)[as.character(clu_genes)%in%as.character(clusters_consistent)]
    p2 = plot_heatmap(list_norm, genes_tmp, clu_genes, pseudotime_rna)
    print(p2)
    dev.off()
  }
 
  
  
  if(interactive){
    print(p1)
    
    print(p2)
  }
  # print(p1)
  # print(p2)
  
}

