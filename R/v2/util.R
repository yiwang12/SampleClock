
# Performs quantile normalization on a matrix of expression values
#
# @param x Matrix of expression values (genes in rows, samples in columns)
# @param preserve_zero Logical indicating whether to preserve zero values (default: TRUE)
#
# @details
# Implementation based on Dave Tang's Blog method:
# 1. Ranks values within each sample
# 2. Takes mean across samples for each rank
# 3. Assigns mean values back based on original ranks
# 4. Optionally preserves zero values from original data 
#
# @return Matrix with same dimensions as input, containing normalized values
#
# @examples
# counts <- matrix(rpois(1000, lambda = 10), ncol = 10)
# normalized_counts <- myQuantileNorm(counts)
#
# @references
# Dave Tang's Blog: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
#
# @export
myQuantileNorm <- function(x, preserve_zero = T){
  # based on Dave Tang's Blog: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
  
  xrank = apply(x,2,rank,ties.method="min")
  xsorted =  apply(x, 2, sort)
  xmean = apply(xsorted, 1, mean)
  
  index_to_mean = function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  xqnorm = apply(xrank, 2, index_to_mean, my_mean= xmean)
  if (preserve_zero){
    xqnorm[which(x == 0)] = 0
  }
  rownames(xqnorm) =  rownames(x)
  return(xqnorm)
}


# Generates pseudobulk data from single-cell expression data by aggregating cells within clusters and samples
#
# Parameters:
# @param s              List of expression matrices or single expression matrix
# @param clu            Vector of cluster labels for each cell
# @param pt             Vector of sample/patient IDs for each cell
# @param topclu         Subset of clusters to analyze (default: all clusters)
# @param batch          Vector of batch labels for batch correction
# @param HVG            Whether to identify highly variable genes (default: TRUE)
# @param qnorm         Whether to perform quantile normalization (default: FALSE)
# @param combat        Whether to perform ComBat batch correction (default: FALSE)
# @param preserve_zero Whether to preserve zeros in quantile normalization (default: TRUE)
#
# Returns:
# List containing:
# - all: Complete pseudobulk expression matrix
# - hvg: Pseudobulk expression matrix for highly variable genes only
#
# @export
get_cluster_pb = function(s, clu, pt,
                          topclu = NULL, batch = NULL, 
                          HVG = T, qnorm = F, combat = F,
                          preserve_zero = T){
  if (is.factor(clu))
    uc = levels(clu)
  else uc = unique(clu) 
  m <- sapply(uc,function(sc) {
    tmp <- s[,which(clu==sc), drop = F]
    p <- pt[which(clu == sc)]
    sapply(unique(p),function(sp) rowSums(tmp[,p==sp,drop=F]))
  })
  d <- lapply(m,function(i) {
    rc <- colSums(i)/1e6
    log2(t(t(i)/rc + 1))
  })
  names(d) = uc
  print(str(d))
  
  ## combat
  if(is.null(topclu))
    topclu = 1:length(s)
  s = s[topclu]
  d <- sapply(1:length(s),function(i) {
    d <- s[[i]][!grepl('MT-',rownames(s[[i]])),]
    d <- d[rowMeans(d > 0.01) > 0.1,]
    
    
    ## qnorm
    if (qnorm){
      dn<- dimnames(d)
      #d<- myQuantileNorm(d, preserve_zero = preserve_zero)
      d = normalize.quantiles(d)
      dimnames(d) = dn; # break ties # RUV? # differential raw pseudo bulk, regress out batch effect?
    }
    
    ## combat
    if (combat){
      match_id = match(colnames(d),names(batch))
      d = ComBat(d,batch = batch[match_id])
    }
    
    ## HVG
    if (HVG){
      cm <- rowMeans(d)
      csd <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
      mod <- loess(csd~cm)
      rid <- which(resid(mod) > 0)
    } else{
      rid = 1:nrow(d)
    }
    
    rownames(d) <- paste0(rownames(d),':',names(s)[i])
        
    return(list(all = d, hvg = d[rid, ]))
  }) 
  
  n <- table(unlist(sapply(d,function(x) colnames(x$all))))
  n <- names(n)[n==length(d)]
  d.all <- do.call(rbind,lapply(d,function(i) i$all[,n]))
  d.hvg <- do.call(rbind,lapply(d,function(i) i$hvg[,n]))
  return(list(all = d.all, hvg = d.hvg))
}

# Generates pseudobulk profiles from single-cell data with preprocessing options
#
# Parameters:
# @param s              Single-cell expression matrix (genes x cells)
# @param pt             Vector of sample/patient IDs for each cell
# @param HVG           Whether to select highly variable genes (default: TRUE)
# @param qnorm         Whether to perform quantile normalization (default: FALSE)
# @param combat        Whether to apply ComBat batch correction (default: FALSE)
# @param batch         Batch labels for ComBat correction
# @param preserve_zero Whether to preserve zeros in quantile normalization (default: TRUE)
#
# Returns:
# List containing:
# - all: Full preprocessed pseudobulk matrix
# - hvg: Subset of highly variable genes (if HVG=TRUE)
#
# @export
get_sample_pb = function(s, pt, HVG = T, qnorm = F, combat = F, batch = NULL,
                         preserve_zero = T){
  m <- sapply(unique(pt),function(sp){
    rowSums(s[,pt ==sp,drop=F])
  })
  rc <- colSums(m)/1e6
  d <- log2(t(t(m)/rc + 1))
  
  ## pre-filtering
  d = d[!grepl('MT-', rownames(d)),]
  d = d[rowMeans(d > 0.01) > 0.1, ]
  
  
  
  ## qnorm
  if (qnorm){
    dn<- dimnames(d)
    #d<- myQuantileNorm(d, preserve_zero = preserve_zero)
    d = normalize.quantiles(d)
    dimnames(d) = dn; # break ties # RUV? # differential raw pseudo bulk, regress out batch effect?
  }
  
  ## combat
  if (combat){
    match_id = match(colnames(d),names(batch))
    d = ComBat(d,batch = batch[match_id])
  }
  
  ## hvg
  if (HVG){
    cm <- rowMeans(d)
    csd <- sqrt((rowMeans(d*d) - cm^2) / (ncol(d) - 1) * ncol(d))
    mod <- loess(csd~cm)
    rid <- which(resid(mod) > 0)
  } else{
    rid = 1:nrow(d)
  }
  return(list(all = d, hvg = d[rid, ]))
}

#' Apply centered log-ratio (CLR) transformation to compositional data
#'
#' @param X Matrix where rows are clusters and columns are samples
#' @param pseudocount Small value added to avoid log(0), defaults to 1
#'
#' @return Matrix of CLR-transformed values
#'
#' @details 
#' The CLR transformation is computed as:
#' CLR(x) = log(x) - mean(log(x))
#' where x is a composition vector
#'
#' @export
clr_transform = function(X,pseudocount = 1){
  # X: each row is a cluster, each column is a sample
  X = X + pseudocount # pseudocount
  ans = apply(X, 2, function(x){log(x/(exp(mean(log(x)))))})
  return(ans)
}


#' Remove batch effects from expression data using linear modeling
#'
#' @param dat Matrix of expression data (features x samples)
#' @param batch Factor indicating batch assignments for each sample
#' @param adj_cov Additional covariates matrix to preserve during batch correction
#'
#' @return Batch-corrected expression matrix with same dimensions as input
#' 
#' @details
#' Implements batch correction by:
#' 1. Fitting linear model with batch effects
#' 2. Estimating batch parameters
#' 3. Subtracting batch effects while preserving biological variation
#'
#' @export
remove_batch_effect = function(dat, batch, adj_cov = NULL){
  batchmod = model.matrix(~ -1 + batch)
  nbatch = ncol(batchmod)
  design = cbind(batchmod, adj_cov)
  
  dat.adj = t(sapply(1:nrow(dat),function(r){
    y = dat[r,]
    beta.hat = solve(crossprod(design), crossprod(design, y))
    gamma.hat = batchmod %*% beta.hat[1:nbatch,]
    grand.mean = mean(gamma.hat)
    gamma.hat = gamma.hat-grand.mean
    y - gamma.hat
  }))
  dimnames(dat.adj) = dimnames(dat)
  return(dat.adj)
}

# Calculate cell type proportions with optional CLR transformation and batch correction
#
# @param clu Vector of cell type labels for each cell
# @param pt Vector of sample/patient IDs corresponding to each cell
# @param clr Logical, whether to apply centered log-ratio transformation (default: FALSE)
# @param batch_correction Logical, whether to perform batch effect correction (default: FALSE)
# @param batch Vector of batch labels for correction
# @param adj_cov Additional covariates matrix for batch correction
#
# @return Matrix of cell type proportions (rows=cell types, columns=samples)
#
# @export
get_ctp = function(clu, pt, clr = F, batch_correction = F, batch = NULL, adj_cov = NULL){
  cell_type_count = table(clu, pt) %>% unclass()
  cell_type_prop = apply(cell_type_count,2,function(x) x/sum(x))
  
  if(clr){
    cell_type_prop = clr_transform(cell_type_count, 1)
  }
  if(batch_correction){
    match_id = match(colnames(cell_type_prop), names(batch))
    cell_type_prop = remove_batch_effect(cell_type_prop, batch[match_id], adj_cov =  NULL)
  }
  return(cell_type_prop)
}

# Evaluates pseudotime ordering quality using different metrics
#
# Parameters:
# @param sample_metadata   Data frame containing sample metadata and pseudotime values
# @param rd               Matrix of reduced dimensions
# @param orth_dist       Vector of orthogonal distances from pseudotime path
# @param sample_name     Column name for sample identifiers (default: 'sample')
# @param pseudotime_name Column name for pseudotime values (default: 'pseudotime')
# @param severity_num_name Column name for numerical severity values (default: 'sev.level')
# @param cluster_name    Column name for cluster assignments (default: 'cluster')
# @param eval_method     Evaluation method to use:
#                       - 'corr': Spearman correlation with severity
#                       - 'pairwise_order': Fraction of correctly ordered sample pairs
#                       - 'R2': Goodness of fit using orthogonal distances
#
# Returns:
# Numeric value representing pseudotime quality:
# - For 'corr': correlation coefficient (-1 to 1)
# - For 'pairwise_order': fraction of correct orderings (0 to 1)
# - For 'R2': explained variance ratio (0 to 1)
#
# @export
get_eval_metric = function(sample_metadata, rd, orth_dist, 
                           sample_name = 'sample',
                           pseudotime_name = 'pseudotime',
                           severity_num_name = 'sev.level',
                           cluster_name = 'cluster',
                           eval_method = c('corr','pairwise_order', 'R2')){
  eval_method = match.arg(eval_method)
  
  if (eval_method == 'corr'){
    ans = cor(sample_metadata[,severity_num_name], sample_metadata[,pseudotime_name],
              use = 'complete')
  }
  
  if (eval_method == 'pairwise_order'){
    pt = sample_metadata[, sample_name]
    submeta = sample_metadata %>% select(!!sample_name, !!pseudotime_name, !! severity_num_name)
    df = expand.grid(x1 = pt, x2 = pt, stringsAsFactors = F) %>%
      filter(x1 != x2)
    df = df %>%
      inner_join(submeta, by = c('x1' = sample_name)) %>%
      inner_join(submeta, by = c('x2' = sample_name)) 
    colnames(df) = c('x1','x2','ptime1','sev1','ptime2','sev2')
    df = df %>%
      mutate(correct_order = ((sev1 <= sev2) == (ptime1 <= ptime2)))
    ans = mean(df$correct_order)
    
  }
  
  if (eval_method == 'R2'){
    tss = sum((sweep(rd, 2, colMeans(rd), '-'))^2)
    rss = sum(orth_dist)
    ans = 1- rss/tss
  }
  return(ans)
}
