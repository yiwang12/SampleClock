
# library(matrixStats)
# 
# ###############################################################
# setwd("/dcl01/hongkai/data/data/ywang/scATAC/cvd/GO_GREAT/")
# 
# ###############################################################
# load("/dcl01/hongkai/data/data/ywang/scATAC/cvd/out_pseudotimeUsingMono/kgamSignificant_means_by_fittedValue_TFmotif/mats_binned_ATAC_RNA.RData")
# # save(mat_rna_binned,mat_atac_binned,file="mats_binned_ATAC_RNA.RData")
# nbin = ncol(mat_rna_binned)
# 
# ###############################################################
# clu_genes=readRDS("clu_genes.rds")
# clu_unique = unique(clu_genes)

#################################### 
#  find gene clusters that are consistent and divergent across cell types

# get_consistent_divergent_clus()


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


# saveRDS(list_consistency_trend_acrossCellTypes,file="list_consistency_trend_acrossCellTypes_gene.rds")

####################### pick top clusters that are consistent

# 
# saveRDS(min_cor_perClu,file="min_cor_perClu_gene.rds")
# 
# 
# saveRDS(clusters_consistent,file="geneClusters_consistent.rds")
# saveRDS(clusters_divergent,file="geneClusters_divergent.rds")

####################### 





