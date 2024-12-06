# library(pheatmap)
# library(ggplot2)
library(matrixStats)
#################################### 
# setwd("/dcl01/hongkai/data/data/ywang/scATAC/cvd/GO_GREAT/")

#################################### functions for data processing
get_processed_peaks <- function(mat_, RNA = T){
  rownames(mat_)=sub(".*[_]","",rownames(mat_))
  if(!RNA){
    c_=colSums(mat_)
    mat_=t(t(mat_)/c_)*10^6
  }

  mat_
}

# mat_atac_aggregated = list_atacPeak_norm[[cellType]]
# pseudotime_atac_agg = pseudotime[colnames(list_atacPeak_norm[[cellType]])]
# list_atacPeak_binned[[cellType]] = get_binned(list_atacPeak_norm[[cellType]], 
                                              # pseudotime[colnames(list_atacPeak_norm[[cellType]])])
# }
get_binned<-function(mat_atac_aggregated, pseudotime_atac_agg){
  pseudotime_atac_agg[is.na(pseudotime_atac_agg)]=max(pseudotime_atac_agg[!is.na(pseudotime_atac_agg)])*2
  
  if(sum(is.na(names(pseudotime_atac_agg)))>0){
    kp_ = which(!is.na(names(pseudotime_atac_agg)))
    pseudotime_atac_agg = pseudotime_atac_agg[kp_]
    mat_atac_aggregated = mat_atac_aggregated[,kp_]
  }

  
  # x <- order(pseudotime_atac_agg, decreasing = F)
  size_bin = round((max(pseudotime_atac_agg)/5))
  
  
  # pseudotime_atac_agg[bins[[1]]]
  
  mat_atac_aggregated_binned=array(dim=c(nrow(mat_atac_aggregated),length(bins)))
  rownames(mat_atac_aggregated_binned) = rownames(mat_atac_aggregated)
  
  for(kbin in 1:length(bins)){
    kst = size_bin*(kbin-1) 
    kend = size_bin*kbin
    if(kbin == length(bins)){kend=max(pseudotime_atac_agg)}
    
    which_withBin = which(pseudotime_atac_agg > kst & pseudotime_atac_agg <= kend)
    
    if(length(which_withBin)>0){
      if(length(which_withBin) > 1){
        mat_atac_aggregated_binned[,kbin] = rowMeans(mat_atac_aggregated[,which_withBin])
      }else{
        mat_atac_aggregated_binned[,kbin] = mat_atac_aggregated[,which_withBin]
      }
      
    }
    
  }
  
 
  rownames(mat_atac_aggregated_binned) = sub(".*[_]", "", (rownames(mat_atac_aggregated_binned)))
  mat_atac_aggregated_binned
}






get_normed_each_cellType <- function(ptime, cellTypes, RNA = T){
  list_atacPeak_norm = list()
  list_atacPeak_binned = list()
  
  # cellTypes = names(ptime$pb.ls$all)
  # cellTypes = cellTypes[cellTypes != "All"]
  for(cellType in cellTypes){
    if(cellType == cellTypes[1]){
      samples_kp = colnames(ptime$pb.ls$all[[cellType]])
    }else{
      samples_kp = intersect(samples_kp, colnames(ptime$pb.ls$all[[cellType]]))
    }
    # samples_kp = c()
    
  }
  
  for(cellType in cellTypes){
    print(cellType)
    # cellType_name = cellType
    # if(cellType=="M"){cellType_name="Mono"}
    list_atacPeak_norm[[cellType]] = get_processed_peaks(ptime$pb.ls$all[[cellType]][,samples_kp], RNA = RNA)
    
    # list_atacPeak_binned[[cellType]] = get_binned(list_atacPeak_norm[[cellType]],
                                                  # pseudotime_atac)
  }
  
  list_atacPeak_norm
}


# list_rnaGene_norm_binned = get_normed_binned_each_cellType(ptime_RNA, cellTypes_, out_$pseudotime_rna, RNA = T)
# ptime = ptime_RNA
# cellTypes = cellTypes_
# pseudotime = out_$pseudotime_rna
# RNA = T

get_normed_binned_each_cellType <- function(ptime, cellTypes, pseudotime, RNA = T){
  
  list_atacPeak_norm = list()
  list_atacPeak_binned = list()
  
  
  for(cellType in cellTypes){
    if(cellType == cellTypes[1]){
      samples_kp = colnames(ptime$pb.ls$all[[cellType]])
    }else{
      samples_kp = intersect(samples_kp, colnames(ptime$pb.ls$all[[cellType]]))
    }
  }
  
  for(cellType in cellTypes){
    print(cellType)
    list_atacPeak_norm[[cellType]] = get_processed_peaks(ptime$pb.ls$all[[cellType]][, samples_kp], RNA = RNA)
    
    list_atacPeak_binned[[cellType]] = get_binned(list_atacPeak_norm[[cellType]], 
                                                  pseudotime[colnames(list_atacPeak_norm[[cellType]])])
  }
  
  list_atacPeak_binned
}





FindHVG_clusters <- function(Zscore_concatenated, list_atacPeak_norm, list_atacPeak_binned, cellTypes, nfeature=10000){
  
  list_Zscore = list()
  for(cellType in cellTypes){
    Zscore_tmp = Zscore_concatenated[grep(paste0(cellType,"."), names(Zscore_concatenated))]
    names(Zscore_tmp) = sub(paste0(cellType, "."), "", names(Zscore_tmp))
    
    list_Zscore[[cellType]] = Zscore_tmp
  }
  
  list_genes_perCellType = list()
  for(cellType in cellTypes){
    list_genes_perCellType[[cellType]] = names(list_Zscore[[cellType]])
  }
  
  #### cluster genes
  ## step1 select HVGs in each cell type
  HVGs_gene = c()
  for(cellType in cellTypes){
    Zscore_tmp = list_Zscore[[cellType]]
    HVGs_gene_tmp  = names(Zscore_tmp)[order(Zscore_tmp, decreasing = T)][1:nfeature]
    # HVGs_gene_tmp = sub(paste0(cellType, "."), "", HVGs_gene_tmp)
    HVGs_gene = unique(c(HVGs_gene,
                         HVGs_gene_tmp))
    
  }
  
  genes_kp = Reduce(intersect,list_genes_perCellType)
  
  # HVGs_gene = Reduce(intersect,genes_kp)
  HVGs_gene = intersect(HVGs_gene,genes_kp)
  
  
  ## step2 cluster genes
  # concatenate the gene exp matrix across cell types
  for(cellType in cellTypes){
    mat_gene_tmp = list_atacPeak_binned[[cellType]]
    # rownames(mat_gene_tmp) = sub(paste0(cellType, "."), "", rownames(mat_gene_tmp))
    
    if(cellType==cellTypes[1]){mat_gene_con = t(scale(t(mat_gene_tmp[HVGs_gene,])))}else{
      mat_gene_con = cbind(mat_gene_con, t(scale(t(mat_gene_tmp[HVGs_gene,]))))
    }
    
  }
  
  set.seed(111)
  clu_genes = kmeans(mat_gene_con, 20)$cluster
  names(clu_genes) = HVGs_gene

  clu_genes
  
}


######


