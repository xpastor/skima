% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differentialMethylation.R
\name{report_significance}
\alias{report_significance}
\title{Report significant results}
\usage{
report_significance(mat, fit, groups = NULL, coef = colnames(fit$design)[2],
  prefix = coef, invisible = T)
}
\arguments{
\item{mat}{matrix used to report the results}

\item{fit}{fitLM object from limma}

\item{groups}{a named vector with the classification used to report results}

\item{coef}{a character with the tested condition to be reported}

\item{prefix}{character, the file prefix to report the results}

\item{invisible}{a boolean to determine if the output of 'missMethyl::gometh' is returned}
}
\value{
limma::topTable
}
\description{
Significant results reported as PCA plots and dendrograms
}
\details{
Extracts the significance analysis from a given variable. Features with an adjusted P-value lower or equal to 0.05 are taken as significant and used to summarize the results with PCA plots and hierarchical clustering.
}
