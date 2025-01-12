% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trendConsistency_geneCluster.R
\name{get_consistent_divergent_clus}
\alias{get_consistent_divergent_clus}
\title{Identify consistent and divergent clusters across cell types based on binned RNA expression}
\usage{
get_consistent_divergent_clus(list_rnaGene_norm_binned, clu_genes, celltypes)
}
\arguments{
\item{list_rnaGene_norm_binned}{List of matrices containing binned normalized RNA expression per cell type}

\item{clu_genes}{Vector of cluster assignments for genes}

\item{celltypes}{Vector of cell type names}
}
\value{
List containing:
\item{clusters_consistent}{Vector of cluster IDs showing consistent patterns}
\item{clusters_divergent}{Vector of cluster IDs showing divergent patterns}
}
\description{
Identify consistent and divergent clusters across cell types based on binned RNA expression
}
\details{
For each cluster:
\enumerate{
\item Calculates average binned expression profiles across genes
\item Computes correlation between cell types
\item Identifies clusters as consistent (min correlation > 0) or divergent (min correlation < 0)
}
}
