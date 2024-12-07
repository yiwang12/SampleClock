#' Process peak matrix data with options for RNA and non-RNA normalization
#'
#' @description
#' Processes a matrix of peak data by cleaning rownames and optionally performing
#' counts per million (CPM) normalization for non-RNA data
#'
#' @param mat_ Input matrix to process
#' @param RNA Logical indicating if data is RNA. If FALSE, performs CPM normalization (default: TRUE)
#'
#' @return Processed matrix with cleaned rownames and optional CPM normalization
#'
#' @details
#' For RNA data, only cleans rownames by removing prefixes
#' For non-RNA data, additionally performs CPM normalization:
#' 1. Calculates column sums
#' 2. Divides each column by its sum
#' 3. Multiplies by 10^6 for CPM
#'
#' @export
get_processed_peaks <- function(mat_, RNA = T){
  rownames(mat_)=sub(".*[_]","",rownames(mat_))
  if(!RNA){
    c_=colSums(mat_)
    mat_=t(t(mat_)/c_)*10^6
  }

  mat_
}
#' Process peak matrix data with options for RNA and non-RNA normalization
#'
#' @param mat_ Input matrix to process
#' @param RNA Logical indicating if data is RNA (default: TRUE)
#'
#' @return Processed matrix with:
#' - Cleaned rownames (prefixes removed)
#' - Optional CPM normalization for non-RNA data
#'
#' @export
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





#' Normalize data for each cell type
#'
#' @param ptime List containing pseudotime analysis results with slot pb.ls$all containing data per cell type
#' @param cellTypes Vector of cell type names to process
#' @param RNA Logical indicating if processing RNA data (TRUE) or ATAC data (FALSE)
#'
#' @return A list where each element contains normalized data for a cell type
#'
#' @details 
#' Function identifies common samples across cell types and normalizes data for each cell type.
#' For RNA data, performs RNA-specific normalizations. For ATAC data, performs ATAC-specific processing.
#'
#' @export
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

#' Process and bin data for each cell type
#'
#' @param ptime List containing pseudobulk data in ptime$pb.ls$all
#' @param cellTypes Vector of cell type names to process
#' @param pseudotime Named vector of pseudotime values for each sample
#' @param RNA Logical indicating if input is RNA data (default: TRUE)
#'
#' @return List of binned data for each cell type after normalization
#'
#' @details
#' 1. Finds common samples across all cell types
#' 2. For each cell type:
#'    - Normalizes data using get_processed_peaks()
#'    - Bins normalized data using get_binned()
#'
#' @export
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



#' Find and cluster highly variable genes across cell types
#'
#' @param Zscore_concatenated Named vector of Z-scores for gene variability, concatenated across cell types
#' @param list_atacPeak_norm List of normalized ATAC peak matrices for each cell type
#' @param list_atacPeak_binned List of binned ATAC peak matrices for each cell type
#' @param cellTypes Vector of cell type names to analyze
#' @param nfeature Number of top variable genes to select per cell type (default: 10000)
#'
#' @return Named vector of cluster assignments for highly variable genes
#'
#' @details
#' The function:
#' 1. Separates Z-scores by cell type
#' 2. Selects top variable genes per cell type
#' 3. Finds genes present across all cell types
#' 4. Creates scaled expression matrix across cell types
#' 5. Performs k-means clustering (k=20) on combined expression
#'
#' @export
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





