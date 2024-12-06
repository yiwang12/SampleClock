library(Matrix)
library(preprocessCore)
library(sva)
library(vegan)
library(cluster)

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


clr_transform = function(X,pseudocount = 1){
  # X: each row is a cluster, each column is a sample
  X = X + pseudocount # pseudocount
  ans = apply(X, 2, function(x){log(x/(exp(mean(log(x)))))})
  return(ans)
}

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
