% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{preprocess}
\alias{preprocess}
\title{Array correction, filtering and QC.}
\usage{
preprocess(idat_dir, targets, blacklist = NULL, batch.vars = NULL,
  removeSNPs = F, population = NULL, usePredictedSex = F, runCNV = F,
  outdir = getwd())
}
\arguments{
\item{idat_dir}{character. Path to the directory with the idat files.}

\item{targets}{data frame with sample annotation}

\item{blacklist}{character. Path to a file with probes to be flagged or removed from the analysis.}

\item{batch.vars}{character. String with a list of columns from \code{'sample.annotation'} separated by comas to be used as confounding variables.}

\item{removeSNPs}{boolean. Should polymorphic probes with an allellic frequency higher than 1/nrow(sample.annotation) in the population given in 'population' be removed?}

\item{population}{character. A string defining the population where to compare allellic frequencies. One of: NULL, '', 'AFR', 'AMR', 'ASN', 'EUR', 'SAS', 'EAS'. NULL and '' will compare to global allellic frequencies.}

\item{usePredictedSex}{boolean. Should the predicted sex be included in the analysis, and used as confounding variable?}

\item{runCNV}{boolean. Should the CNV analysis be run?}

\item{outdir}{character. Directory where all the results will be writen.}
}
\value{
\describe{A list with the following elements:
  \item{M:}{A matrix with M values}
  \item{betas:}{A matrix with Beta values}
  \item{targets:}{A data frame with samples annotation}
  \item{array.annot.gr:}{A GRanges object with array annotation}
}
}
\description{
Function to prepare the data for further analysis, normalizes the arrays, flags cross-reactive and polymorphic probes, removes low quality values and runs QC for raw and normalized data.
}
\details{
'targets' is a data frame that can be derived from Illumina's sample sheet. The columns 'Sample_Name', 'Sentrix_ID' and 'Sentrix_Position' are required. The columns 'Sample_Well', 'Sample_Plate', 'Sample_Group' and 'Pool_ID' are ignored. Any other column will be used as confounding factor or as variable to analyse. The values in 'Sample_Name' must be unique sample identifiers.

'blacklist' contains a list of probes to be flagged or removed from the analysis. Each row may have one or two columns, separated by a tab and the first one being the probe 'cg' identifier. If there is no other column the probe will be flagged. In case there is a second column, it contains the sample identifier, and the value for that probe in that sample will be replaced by 'NA'. This file can be used to indicate polymorphic probes identified by sequencing.

In case polymorphic probes can not be inferred, 'removeEuropeanSNPs' indicates if to flag polymorphic probes based on allellic frequencies in european populations. Probes containing a polymorphism that could be present in at least one of the samples given the allellic frequency in european populations will be flagged.

The gender of the samples is always predicted and reported. 'usePredictedSex' is used to specify if the preditcted gender should be used as confounding variable.
}
\author{
Xavier Pastor \email{x.pastorhostench@dkfz.de}
}
