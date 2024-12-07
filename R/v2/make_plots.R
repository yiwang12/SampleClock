

# Generates visualization plots for integrated RNA and ATAC-seq sample data
#
# Parameters:
# @param out_           List containing integrated analysis results with components:
#                      - integrated: Seurat object with RNA/ATAC integration
#                      - pseudotime_rna: RNA pseudotime values
#                      - pseudotime_atac: ATAC pseudotime values
#                      - CC1, CC2: Canonical correlation components
#                      - slope: Slope for projection line
#                      - pseudotime_: Combined pseudotime values
# @param interactive    Logical, whether to display plots interactively (default: FALSE)
# @param dir_output    Directory path for saving plots. If NULL, plots aren't saved
#
# Details:
# Creates 4 visualization plots:
# 1. PCA colored by severity
# 2. CCA colored by severity 
# 3. PCA colored by pseudotime
# 4. CCA colored by pseudotime
# Each plot shows RNA and ATAC modalities side by side
#
# @import Seurat
# @import ggplot2
# @import grDevices
#
# @export
plot_integrated_samples <- function(out_, interactive = F, dir_output = NULL){
  
  # pseudotime_atac" "pseudotime_rna
  severity.color = c("#4DAF4A" ,"#377EB8",  "orange", "#E41A1C" )
  names(severity.color)[1:4]=c("Healthy","Mild","Moderate","Severe")
  #####
  integrated=out_$integrated
  # head(integrated@meta.data)
  # head(integrated@meta.data)
  integrated@meta.data$severity=as.character(integrated@meta.data$sev[,1])
  integrated@meta.data$modality=c(rep("RNA",length(out_$pseudotime_rna)),
                                  rep("ATAC",length(out_$pseudotime_atac)))
  
  
  pc=Embeddings(integrated[,], reduction = "pca")
  slope=out_$slope
  
  # kp_ = which(-pc[,1]<30)
  # sum(-pc[,1]>=30)
  # 12
  # p = DimPlot(integrated[,kp_], 
  #           group.by = "severity",
  #           # cols=alpha(my_cols,0.01),
  #           split.by="modality",
  #           ncol=2,
  #           label=F)#+ggtitle("cell type")
  # p$layers[[1]]$aes_params$alpha =  .3
  # p
  
  severity_integrated=integrated@meta.data$severity
  
  print("severity_integrated")
  
  print(severity_integrated)
  
  severity_integrated[severity_integrated=="0"]="Healthy"
  severity_integrated[severity_integrated=="1"]="Mild"
  severity_integrated[severity_integrated=="2"]="Moderate"
  severity_integrated[severity_integrated=="3"]="Severe"
  severity_integrated=factor(severity_integrated,levels=c("Healthy","Mild","Moderate","Severe"))
  print(severity_integrated)
  
  sample_info=data.frame(PC1=pc[,1],
                         PC2=pc[,2],
                         CC1=out_$CC1,
                         CC2=out_$CC2,
                         severity=severity_integrated,
                         slope=slope,
                         nFeature_RNA=integrated@meta.data$nFeature_RNA,
                         pseudotime=out_$pseudotime_,
                         modality= integrated@meta.data$modality)
  
  
  # saveRDS(sample_info,file= paste0(dir_output,
            # "/sample_info_sampleLevelIntegrated_MonoOnly.rds"))
  
  # saveRDS(slope,file=
            # "/dcl01/Chongkai/data/data/ywang/scATAC/cvd/out_pseudotimeUsingMono/slope_sampleLevelIntegrated_MonoOnly.rds")
  
  
  ##### create plots
  # sample_info = readRDS("/dcl01/hongkai/data/data/ywang/scATAC/cvd/out_pseudotimeUsingMono/sample_info_sampleLevelIntegrated_MonoOnly.rds")
  # slope = readRDS("/dcl01/hongkai/data/data/ywang/scATAC/cvd/out_pseudotimeUsingMono/slope_sampleLevelIntegrated_MonoOnly.rds")
  
  # kp_ = which(sample_info$PC1> -30)
  # sample_info=sample_info[kp_,]
  
  #####
  if(!is.null(plot_integrated_samples)){
    dir.create(dir_output)
    
    pdf(paste0(dir_output, "/SampleLevelIntegration_MonoOnly_PCA_colBySevs.pdf"),width=7,height=4)
    p1 <- ggplot(aes(x = PC1, y = PC2, color = severity), data= sample_info)+
      geom_point(size=.3)+
      scale_colour_manual(values = severity.color)+
      geom_abline(slope = slope, intercept = 0, linetype = 'dashed', col = 'black')+
      labs(color = "Severity") +
      ggtitle("Severity")+ 
      facet_wrap(~modality,ncol=2)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(p1)
    dev.off()
    
    pdf(paste0(dir_output, "/SampleLevelIntegration_MonoOnly_CCA_colBySev.pdf"),width=7,height=4)
    p2 <- ggplot(aes(x = CC1, y = CC2, color = severity), data= sample_info)+
      geom_point(size=.3)+
      scale_colour_manual(values = severity.color)+
      labs(color = "Severity") +
      ggtitle("Severity")+ 
      facet_wrap(~modality,ncol=2)+
      
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(p2)
    dev.off()
    
    pdf(paste0(dir_output, "/SampleLevelIntegration_MonoOnly_PCA_colByPseudotime.pdf"),width=7,height=4)
    p3 <- ggplot(aes(x = PC1, y = PC2, color = pseudotime), data= sample_info)+
      geom_point(size=.3)+
      geom_abline(slope = slope, intercept = 0, linetype = 'dashed', col = 'black')+
      ggtitle("Pseudotime")+   facet_wrap(~modality,ncol=2)+
      
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(p3)
    dev.off()
    
    pdf(paste0(dir_output, "/SampleLevelIntegration_MonoOnly_CCA_colByPseudotime.pdf"),width=7,height=4)
    p4 <-  ggplot(aes(x = CC1, y = CC2, color = pseudotime), data= sample_info)+
      geom_point(size=.3)+
      ggtitle("Pseudotime")+   facet_wrap(~modality,ncol=2)+
      
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(p4)
    dev.off()
  }
  
  if(interactive){
    print(p1)
    
    print(p3)
  }
  # print(p1)
  # print(p2)
  # print(p3)
  # print(p4)
  
}



