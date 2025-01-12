% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_genes_heatmap.R
\name{plot_consist_div}
\alias{plot_consist_div}
\title{Function to generate and save heatmap visualizations for divergent and consistent gene clusters}
\usage{
plot_consist_div(
  list_norm,
  clu_genes,
  clusters_consistent,
  clusters_divergent,
  pseudotime_rna,
  dir_plot = NULL,
  interactive = F
)
}
\arguments{
\item{list_norm}{Normalized expression data list}

\item{clu_genes}{Named vector of gene cluster assignments}

\item{clusters_consistent}{Vector of consistent cluster IDs}

\item{clusters_divergent}{Vector of divergent cluster IDs}

\item{pseudotime_rna}{Pseudotime ordering for samples}

\item{dir_plot}{Output directory for saved plots (default: NULL)}

\item{interactive}{Whether to display plots interactively (default: FALSE)

Returns:
Creates PDF files of heatmaps and optionally displays them interactively}
}
\description{
Parameters:
}
