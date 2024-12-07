

#' Identify consistent and divergent clusters across cell types based on binned RNA expression
#'
#' @param list_rnaGene_norm_binned List of matrices containing binned normalized RNA expression per cell type
#' @param clu_genes Vector of cluster assignments for genes
#' @param celltypes Vector of cell type names
#'
#' @return List containing:
#'   \item{clusters_consistent}{Vector of cluster IDs showing consistent patterns}
#'   \item{clusters_divergent}{Vector of cluster IDs showing divergent patterns}
#'
#' @details 
#' For each cluster:
#' 1. Calculates average binned expression profiles across genes
#' 2. Computes correlation between cell types
#' 3. Identifies clusters as consistent (min correlation > 0) or divergent (min correlation < 0)
#'
#' @export
get_consistent_divergent_clus <- function(list_rnaGene_norm_binned, clu_genes, celltypes){
  
  clu_unique = unique(clu_genes)
  
  list_consistency_trend_acrossCellTypes = list()
  
  nbin = ncol(list_rnaGene_norm_binned[[1]])
  
  for(kclu in 1:length(clu_unique)){
    
    # print(kclu)
    clu_tmp = clu_unique[kclu]
    
    genes_cluster_tmp = names(clu_genes)[clu_genes==clu_tmp]
    
    mat_ave_bineed = array(dim=c(nbin, length(celltypes)))
    colnames(mat_ave_bineed) = celltypes
    
    
    for(celltype in celltypes){
      mat_rna_binned_tmp = list_rnaGene_norm_binned[[celltype]][genes_cluster_tmp, ]
      
      # mat_rna_binned_tmp = list_rnaGene_norm[[celltype]][paste0(celltype,".",genes_cluster_tmp),]
      mat_ave_bineed[, celltype] = colSums(mat_rna_binned_tmp)
    }
    
    list_consistency_trend_acrossCellTypes[[as.character(clu_tmp)]] = cor(mat_ave_bineed)
    
  }
  
  
  min_cor_perClu = array(dim=length(clu_unique))
  names(min_cor_perClu) = as.character(clu_unique)
  for(clu_tmp in as.character(clu_unique)){
    cor_tmp = list_consistency_trend_acrossCellTypes[[as.character(clu_tmp)]]
    min_cor_perClu[clu_tmp] = min(cor_tmp[upper.tri(cor_tmp)])
  }
  
  # sort(min_cor_perClu)
  clusters_consistent = names(min_cor_perClu)[min_cor_perClu>0]
  clusters_divergent = names(min_cor_perClu)[min_cor_perClu<0]
  
  list(clusters_consistent = clusters_consistent,
       clusters_divergent = clusters_divergent)
}








