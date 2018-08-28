#' Execution of SKIMA.
#'
#' Runs the SKIMA workflow.
#'
#' @param idat_dir character. Path to the directory with the idat files.
#' @param sample.annotation character. Path to the sample sheet with all the annotation about the samples. Must have a column named 'Sample_Name' with unique sample identifiers, a column named 'Sentrix_ID' with the array identifier and a column named 'Sentrix_Position' with the position in the array where the sample was hybridized.
#' @param blacklist character. Path to a file with probes to be flagged or removed from the analysis.
#' @param outdir character. Directory where all the results will be writen.
#' @param runCNV boolean. Should the CNV analysis be run?
#' @param diffMeth boolean. Should the differential methylation analsysis be run?
#' @param removeSNPs boolean. Should polymorphic probes with an allellic frequency higher than 1/nrow(sample.annotation) in the population given in 'population' be removed?
#' @param population character. A string defining the population where to compare allellic frequencies. One of: NULL, '', 'AFR', 'AMR', 'ASN', 'EUR', 'SAS', 'EAS'. NULL and '' will compare to global allellic frequencies.
#' @param usePredictedSex boolean. Should the predicted sex be included in the analysis, and used as confounding variable?
#' @param batch.vars character. String with a list of columns from \code{'sample.annotation'} separated by comas to be used as confounding variables.
#' @param seed integer. Seed to avoid randomness and make the analysis reproducible.
#' @param ncores integer. Number of cores to use for the multithreaded processes.
#'
#' @return None
#'
#' @author Xavier Pastor \email{x.pastorhostench@@dkfz.de}
#'
#' @examples
#' run()
#'
#' @export
#'
#' @importFrom minfi read.metharray.sheet
run <- function(
	idat_dir = Sys.getenv("IDAT_DIR"),
	sample.annotation = Sys.getenv("SAMPLE_ANNOTATION"),
	blacklist = Sys.getenv("BLACKLIST"),
	outdir = Sys.getenv("OUTDIR"),
	runCNV = as.logical(Sys.getenv("RUN_CNV")),
	diffMeth = as.logical(Sys.getenv("RUN_DIFFERENTIAL_METHYLATION")),
	removeSNPs = as.logical(Sys.getenv("REMOVE_SNPS")),
	population = Sys.getenv("POPULATION"),
	usePredictedSex = as.logical(Sys.getenv("USE_PREDICTED_SEX")),
	batch.vars = Sys.getenv("BATCH_VARS"),
	seed = as.integer(Sys.getenv("SEED")),
	ncores = as.integer(Sys.getenv("NCORES")))
{
	idat_dir <- .empty2null(idat_dir)
	sample.annotation <- .empty2null(sample.annotation)
	blacklist <- .empty2null(blacklist)
	outdir <- .empty2null(outdir)
#	seed <- as.integer(seed)
#	ncores <- as.integer(ncores)
	batch.vars <- .empty2null(batch.vars)
	batch.vars <- if (!is.null(batch.vars)) {if (length(batch.vars) == 1) unlist(strsplit(batch.vars, ','))}
	if (is.null(population)) population <- ''

	# Set default values for empty inputs
	if (is.null(outdir)) outdir <- getwd()
	if (is.na(runCNV)) runCNV <- FALSE
	if (is.na(diffMeth)) diffMeth <- FALSE
	if (is.na(removeSNPs)) removeSNPs <- FALSE
	if (is.na(usePredictedSex)) usePredictedSex <- TRUE
	if (is.na(seed)) seed <- 13
	## Adjust the number of usable cores
	if (Sys.getenv('PBS_NUM_PPN') != '') {
		ncores <- min(ncores, as.integer(Sys.getenv('PBS_NUM_PPN')))
	} else {
		ncores <- min(ncores, parallel::detectCores()-1, na.rm=T)
	}
	ncores <- max(1, ncores, na.rm=T)

	.run_checks()

	### Preparing and checking targets data frame ###
	targets <- if(is.null(sample.annotation)) read.metharray.sheet(idat_dir) else .read.targets(idat_dir, sample.annotation)
	if (!all(batch.vars %in% colnames(targets))) {
		stop("\n\tOne or more of your batch variables are not present in the sample annotation file.")
	}
	if (diffMeth) {
		interest.vars <- .get_interest.vars(targets, batch.vars, FALSE)
		no.replicates <- .check_targets(targets[,interest.vars])
	}
	
	# Define the seed
	set.seed(seed)
	#o#

	methods.txt <- file.path(outdir, 'methods.txt')
	citations.txt <- file.path(outdir, 'citations.txt')
	
	library(tools)
	obj <- preprocess(idat_dir = idat_dir, targets = targets, blacklist = blacklist, batch.vars = batch.vars, removeSNPs = removeSNPs, population = population, usePredictedSex = usePredictedSex, runCNV = runCNV, outdir = outdir)

	M <- obj$M
	betas <- obj$betas
	targets <- obj$targets
	array.annot.gr <- obj$array.annot.gr

#	if (runCNV) CNV()
	if (diffMeth) {
		interest.vars <- .get_interest.vars(targets, batch.vars, usePredictedSex)
		if (isEmpty(interest.vars)) {
			message("There's no factor eligible for a differential methylation analysis and it will be skipped.")
		} else {
			ok <- names(array.annot.gr)[array.annot.gr$score==0]
			fit <- differentialMethylation(
					M=M[ok,],
					targets=targets[,unique(c(batch.vars, interest.vars)),drop=F],
					array.annot.gr=array.annot.gr,
					betas=betas[ok,],
					interest.vars=interest.vars,
					outdir=outdir,
					ncores=ncores,
					invisible=F)
		}
	}
	array.type <- guessArrayType(rownames(betas))
	if (population == '') population <- 'mixed'
	af <- NULL
	if(removeSNPs) af <- 1/ncol(M)
	write.methods(array.type, blacklist, af, population, diffMeth=diffMeth & (length(interest.vars) > 0), outdir)
	write.citations(array.type, diffMeth, outdir)
	report_session(outdir)
}

.empty2null <- function(x)
{
	if (!is.null(x)) {
		if (x == '') x <- NULL
	}
	return(x)
}

.read.targets <- function(idat_dir, sample.annotation)
{
	targets <- read.metharray.sheet(dirname(sample.annotation), paste0('^',basename(sample.annotation),'$'), recursive=F)
	targets$Basename <- paste(targets$Slide, targets$Array, sep='_')
	targets$Basename <- file.path(idat_dir, targets$Basename)
	return(targets)
}

.get_interest.vars <- function(targets, batch.vars, usePredictedSex=F)
{
	# Extract variables of interest from sample sheet
	interest.vars <- colnames(targets)
	noisy.vars <- c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Pool_ID', 'Slide', 'Array', 'Basename', 'filenames')
	interest.vars <- interest.vars[! interest.vars %in% c(noisy.vars, batch.vars)]
	if (usePredictedSex) {
		interest.vars <- c(interest.vars, 'predictedSex')
	}
	return(interest.vars)
}
