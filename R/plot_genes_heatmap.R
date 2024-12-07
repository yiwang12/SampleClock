





#' Create pseudotime-ordered gene expression heatmaps
#'
#' @param list_norm List of normalized expression matrices for different cell types
#' @param genes_ Vector of gene names to plot
#' @param clu_genes_ Named vector of cluster assignments for genes
#' @param pseudotime_ Named vector of pseudotime values for samples
#'
#' @details
#' The function:
#' - Orders samples by pseudotime
#' - Orders genes by cluster
#' - Scales expression values to z-scores
#' - Caps z-scores at [-1,1]
#' - Creates separate heatmaps for each cell type
#' - Annotates genes by cluster assignment
#'
#' @importFrom pheatmap pheatmap
#'
#' @export
plot_heatmap<-function(list_norm, genes_, clu_genes_, pseudotime_){
  my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
               "#E16462FF", "#FCA636FF", "#F0F921FF")
  
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


#' Function to generate and save heatmap visualizations for divergent and consistent gene clusters
#'
#' Parameters:
#' @param list_norm            Normalized expression data list
#' @param clu_genes           Named vector of gene cluster assignments
#' @param clusters_divergent   Vector of divergent cluster IDs
#' @param clusters_consistent  Vector of consistent cluster IDs
#' @param pseudotime_rna      Pseudotime ordering for samples
#' @param dir_plot           Output directory for saved plots (default: NULL)
#' @param interactive        Whether to display plots interactively (default: FALSE)
#'
#' Returns:
#' Creates PDF files of heatmaps and optionally displays them interactively
#'
#' @export
plot_consist_div <- function(list_norm, clu_genes, clusters_consistent, clusters_divergent, 
                             pseudotime_rna, dir_plot=NULL, interactive=F) {
  # Input validation
  if (is.null(list_norm) || length(list_norm) == 0) {
    stop("list_norm cannot be empty or NULL")
  }
  if (is.null(clu_genes) || length(clu_genes) == 0) {
    stop("clu_genes cannot be empty or NULL")
  }
  
  # Create directory if it doesn't exist and dir_plot is provided
  if (!is.null(dir_plot)) {
    dir.create(dir_plot, showWarnings = FALSE)
  }
  
  # Define colors
  my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
               "#E16462FF", "#FCA636FF", "#F0F921FF")
  
  # Initialize plot variables
  p1 <- p2 <- NULL
  
  # Create divergent genes plot
  if (!is.null(dir_plot)) {
    genes_div <- names(clu_genes)[as.character(clu_genes) %in% as.character(clusters_divergent)]
    if (length(genes_div) > 0) {
      pdf(file.path(dir_plot, "geneCluster_heatmap_div.pdf"))
      p1 <- plot_heatmap(list_norm, genes_div, clu_genes, pseudotime_rna)
      print(p1)
      dev.off()
    }
    
    # Create consistent genes plot
    genes_cons <- names(clu_genes)[as.character(clu_genes) %in% as.character(clusters_consistent)]
    if (length(genes_cons) > 0) {
      pdf(file.path(dir_plot, "geneCluster_heatmap_consistent.pdf"))
      p2 <- plot_heatmap(list_norm, genes_cons, clu_genes, pseudotime_rna)
      print(p2)
      dev.off()
    }
  }
  
  # Display plots if interactive is TRUE and plots were created
  if (interactive) {
    if (!is.null(p1)) print(p1)
    if (!is.null(p2)) print(p2)
  }
  
  # Return the plot objects invisibly
  invisible(list(divergent = p1, consistent = p2))
}

# # Helper function to check if input data is valid
# check_plot_data <- function(x) {
#   if (is.null(x) || length(x) == 0 || all(is.na(x))) {
#     return(FALSE)
#   }
#   return(TRUE)
# }
# 
# plot_consist_div <- function(list_norm, clu_genes, clusters_divergent, clusters_consistent, pseudotime_rna, dir_plot=NULL, interactive=F){
#   dir.create(dir_plot)
#   my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
#                "#E16462FF", "#FCA636FF", "#F0F921FF")
#   
#   if(!is.null(dir_plot)){
#     pdf(paste0(dir_plot, "/geneCluster_heatmap_div.pdf"))
#     genes_tmp = names(clu_genes)[as.character(clu_genes)%in%as.character(clusters_divergent)]
#     p1 = plot_heatmap(list_norm, genes_tmp, clu_genes, pseudotime_rna)
#     print(p1)
#     dev.off()
#     
#     
#     pdf(paste0(dir_plot, "/geneCluster_heatmap_consistent.pdf"))
#     genes_tmp = names(clu_genes)[as.character(clu_genes)%in%as.character(clusters_consistent)]
#     p2 = plot_heatmap(list_norm, genes_tmp, clu_genes, pseudotime_rna)
#     print(p2)
#     dev.off()
#   }
#  
#   
#   
#   if(interactive){
#     print(p1)
#     
#     print(p2)
#   }
#   # print(p1)
#   # print(p2)
#   
# }
# 
