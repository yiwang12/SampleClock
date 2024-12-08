% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_overwrite_ptime_1.R
\name{run_pseudo_diff_gene}
\alias{run_pseudo_diff_gene}
\title{Identify differentially expressed genes along pseudotime trajectory}
\usage{
run_pseudo_diff_gene(
  x,
  sample_metadata = NULL,
  cluster = "all",
  sample_name = "sample",
  pseudotime_name = "pseudotime",
  gene_mapping = F,
  deg_parallel = F,
  deg_num_cores = 10
)
}
\arguments{
\item{x}{Matrix of expression values (genes x samples)}

\item{sample_metadata}{Data frame with sample metadata including pseudotime values}

\item{cluster}{Cell/sample cluster identifier (default: 'all')}

\item{sample_name}{Column name for sample IDs in metadata (default: 'sample')}

\item{pseudotime_name}{Column with pseudotime values in metadata (default: 'pseudotime')}

\item{gene_mapping}{Whether to map gene symbols to IDs using org.Hs.eg.db (default: FALSE)}

\item{deg_parallel}{Whether to use parallel processing (default: FALSE)}

\item{deg_num_cores}{Number of cores for parallel processing (default: 10)}
}
\value{
Data frame containing:
\itemize{
\item gene: gene identifiers
\item pval: raw p-values from GAM model fit
\item fdr: adjusted p-values (FDR)
\item effect_size: magnitude of expression changes
\item cell_cluster: cluster identifier
\item signif: boolean indicating significance (FDR < 0.05)
If gene_mapping=TRUE, additional columns:
\item ENSEMBL: Ensembl gene IDs
\item ENTREZID: Entrez gene IDs
}
}
\description{
Identify differentially expressed genes along pseudotime trajectory
}
\details{
Uses generalized additive models (GAM) to test for gene expression changes
along pseudotime. Expression patterns are modeled using cubic splines.
Effect size is calculated as max-min fitted values normalized by residual std dev.
}
