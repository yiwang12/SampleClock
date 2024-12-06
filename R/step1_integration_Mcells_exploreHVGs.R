library(RColorBrewer)
library(ggrepel)
library(Matrix)
library(parallel)
library(igraph)
library(Seurat)
library(ggplot2)


###############################################################
## process data for each cell type 
getSeuratProcessed<-function(mat_log2,nfeature_){
  obj_=CreateSeuratObject(mat_log2)
  obj_=NormalizeData(obj_)
  obj_$RNA@data=mat_log2
  obj_=FindVariableFeatures(obj_,nfeature=nfeature_)
  
  obj_
}


getZscore<-function(obj_){
  z_score=HVFInfo(obj_)$variance.standardized
  names(z_score)=rownames(obj_)
  z_score
}

get_consistency_pseudotime_sev_integratedData<-function(Features_sele_,
                                                        mca.RNA_,
                                                        mca.ATAC_,
                                                        sev_num_rna_,
                                                        sev_num_atac_){
  
  sev_num_rna_ = as.integer(sev_num_rna_)
  sev_num_atac_ = as.integer(sev_num_atac_)
  
  ###
  mca.RNA_ = mca.RNA_[, !is.na(sev_num_rna_)]
  mca.ATAC_ = mca.ATAC_[, !is.na(sev_num_atac_)]
  
  sev_num_rna_ = sev_num_rna_[!is.na(sev_num_rna_)]
  sev_num_atac_ = sev_num_atac_[!is.na(sev_num_atac_)]
  
  ###
  
  
  pancreas.anchors <- FindIntegrationAnchors(object.list = list(RNA=mca.RNA_,
                                                                ATAC=mca.ATAC_),
                                             # reference=2,
                                             # reduction="rpca",
                                             dims=1:5,
                                             k.anchor = 5,
                                             k.filter = 10,
                                             k.score = 10)
  

  pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, 
                                       dims = 1:5, k.weight=5,
                                       features.to.integrate=Features_sele_)
  print(4)
  
  DefaultAssay(pancreas.integrated) <- "integrated"
  print(5)
  
  # Run the standard workflow for visualization and clustering
  pancreas.integrated@assays$integrated@data@x[is.na(pancreas.integrated@assays$integrated@data@x)] <- 0
  pancreas.integrated@assays$integrated@data[is.na(pancreas.integrated@assays$integrated@data)] <- 0
  

  # pancreas.integrated <- NormalizeData(object = pancreas.integrated)
  # pancreas.integrated <- FindVariableFeatures(object = pancreas.integrated)
  VariableFeatures( pancreas.integrated) = Features_sele_
  pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
  print(6)
  
  
  pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 5, verbose = T,
                                features=Features_sele_)
  
  ################## 
  ################## load meta data
  ################## 
  
  print(7)
  
  sev_level=data.frame(sev=c(sev_num_rna_,sev_num_atac_))
  
  rownames(sev_level)=colnames(pancreas.integrated)
  pancreas.integrated=AddMetaData(pancreas.integrated,sev_level)
  
  pancreas.integrated@meta.data$sev=sev_level
  
  print(8)
  print(sev_level)
  print(88)
  
  ##################
  ### CC1
  print(head(pc))
  
  pc=Embeddings(pancreas.integrated[,], reduction = "pca")
  
  severity_num =pancreas.integrated@meta.data$sev[]
  res = cancor(pc,severity_num)
  
  
  
  print(9)
  
  proj = sweep(pc, 2, res$xcenter, '-') %*% res$xcoef[,1, drop = F] #?? CC1
  
  CC1 = proj
  CC2 = sweep(pc, 2, res$xcenter, '-') %*% res$xcoef[,2, drop = F] #?? CC1
  
  slope = res$xcoef["PC_2",1]/res$xcoef["PC_1",1] #???
  print(10)
  
  pseudotime_ = rank(proj)
  pseudotime_=as.integer(pseudotime_)
  print(11)
  
  if (cor(pseudotime_, severity_num)<0)
    pseudotime_ = length(pseudotime_)+1 -pseudotime_
  
  names(pseudotime_)=rownames(CC1)
  
  pseudotime_rna=pseudotime_[colnames(mca.RNA_)]
  pseudotime_atac=pseudotime_[colnames(mca.ATAC_)]
  
  # out_inte$integrated@meta.data$sev=as.character((out_inte$integrated@meta.data$sev)$sev)
  print(12)
  
  # # print(severity_num)
  # print(pseudotime_)
  print(pseudotime_)
  # cor(pseudotime_,(severity_num))
  return(list(cor_=cor(pseudotime_,(severity_num)),
              integrated=pancreas.integrated,
              pseudotime_rna=pseudotime_rna,
              pseudotime_atac=pseudotime_atac,
              pseudotime_=pseudotime_,
              sev_level=sev_level,
              sev_num_rna=sev_num_rna_,
              sev_num_atac=sev_num_atac_,
              # sev_level=data.frame(sev=c(sev_num_rna_,sev_num_atac_))
              # pancreas.integrated=pancreas.integrated,
              res=res,#xcoef,xcenter
              CC1=CC1,
              CC2=CC2,
              slope=slope
  ))
}

# input: SampleClock output object
# output: seurat object of concatenated normalized matrix (log_cpm)
process_pseudobulk <- function(ptime, cellTypes){
  
  # cellTypes = names(ptime$pb.ls$all)
  # cellTypes = cellTypes[cellTypes!="All"]
  
  for(cellType in cellTypes){
    
    pseudobulk_cellType_tmp = ptime$pb.ls$all[[cellType]]
    rownames(pseudobulk_cellType_tmp) = paste0(cellType, ".", rownames(pseudobulk_cellType_tmp))
    pseudobulk_cellType_tmp_norm = t(t(pseudobulk_cellType_tmp)/colSums(pseudobulk_cellType_tmp))*mean(colSums(pseudobulk_cellType_tmp))
    # pseudobulk_cellType_tmp_norm = t(t(pseudobulk_cellType_tmp)/colSums(pseudobulk_cellType_tmp)) * 10^6
    
    
    if(cellType == cellTypes[1]){ 
      samples_kp = colnames(pseudobulk_cellType_tmp_norm) 
      mat_atac = pseudobulk_cellType_tmp_norm
    }else{
      samples_kp = intersect(samples_kp, colnames(pseudobulk_cellType_tmp_norm))
      mat_atac = rbind(mat_atac[,samples_kp], pseudobulk_cellType_tmp_norm[,samples_kp])
    }
    
  }
  
  obj.ATAC_Con <- getSeuratProcessed(mat_atac[,samples_kp],10000)
  
  obj.ATAC_Con
}

# dim(obj.ATAC_Con)
# whihc_M = grep("Mono.", rownames(obj.ATAC_Con))
# obj.ATAC_Con[whihc_M[473:473],]
# obj.ATAC_Con[whihc_M[473],]
# rownames(obj.ATAC_Con)[whihc_M[473]]
# # Mono.Y_RNA
# obj.ATAC_Con = obj.ATAC_Con[,]

# obj.RNA_Con = obj.RNA_Con[grep("Mono.", rownames(obj.RNA_Con)),]

select_HVG_for_integration <- function(obj.RNA_Con, obj.ATAC_Con){
  z_ATAC=getZscore(obj.ATAC_Con)
  
  z_RNA=getZscore(obj.RNA_Con)
  # z_BMTA_atac = c(z_M_ATAC)
  HVGs_atac = names(z_ATAC)[order(z_ATAC,decreasing = T)][1:(min(length(z_ATAC),10000))]
  HVGs_rna = names(z_RNA)[order(z_RNA,decreasing = T)][1:(min(length(z_RNA),5000))]
  # HVGs = intersect(HVGs_atac,HVGs_rna)
  
  # z_RNA_M = z_RNA[grep("Mono.", names(z_RNA))]
  # HVGs_rna_M = names(z_RNA_M)[order(z_RNA_M,decreasing = T)][1:5000]
  # HVGs = intersect(names(z_ATAC), HVGs_rna_M)
  
  HVGs = intersect(names(z_ATAC), HVGs_rna)
  
  HVGs
}
# z_ATAC_M = z_ATAC[grep("Mono.", names(z_ATAC))]
