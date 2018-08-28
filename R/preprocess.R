#' Array correction, filtering and QC.
#'
#' Function to prepare the data for further analysis, normalizes the arrays, flags cross-reactive and polymorphic probes, removes low quality values and runs QC for raw and normalized data.
#' 
#' @details
#' 'targets' is a data frame that can be derived from Illumina's sample sheet. The columns 'Sample_Name', 'Sentrix_ID' and 'Sentrix_Position' are required. The columns 'Sample_Well', 'Sample_Plate', 'Sample_Group' and 'Pool_ID' are ignored. Any other column will be used as confounding factor or as variable to analyse. The values in 'Sample_Name' must be unique sample identifiers.
#'
#' 'blacklist' contains a list of probes to be flagged or removed from the analysis. Each row may have one or two columns, separated by a tab and the first one being the probe 'cg' identifier. If there is no other column the probe will be flagged. In case there is a second column, it contains the sample identifier, and the value for that probe in that sample will be replaced by 'NA'. This file can be used to indicate polymorphic probes identified by sequencing.
#'
#' In case polymorphic probes can not be inferred, 'removeEuropeanSNPs' indicates if to flag polymorphic probes based on allellic frequencies in european populations. Probes containing a polymorphism that could be present in at least one of the samples given the allellic frequency in european populations will be flagged.
#'
#' The gender of the samples is always predicted and reported. 'usePredictedSex' is used to specify if the preditcted gender should be used as confounding variable.
#'
#' @inheritParams run
#' @param targets data frame with sample annotation
#'
#' @return
#' \describe{A list with the following elements:
#'   \item{M:}{A matrix with M values}
#'   \item{betas:}{A matrix with Beta values}
#'   \item{targets:}{A data frame with samples annotation}
#'   \item{array.annot.gr:}{A GRanges object with array annotation}
#' }
#'
#' @author Xavier Pastor \email{x.pastorhostench@@dkfz.de}
#'
#' @export
preprocess <- function(idat_dir, targets, blacklist=NULL, batch.vars=NULL, removeSNPs=F, population=NULL, usePredictedSex=F, runCNV=F, outdir=getwd())
{
	# Load libraries
	library(minfi)
	library(ENmix)
	library(GenomicRanges)

	#### Reading in data ####
	### Reading methylation data ###
	message("Reading in methylation files...")
	rgset <- read.metharray.exp(targets=targets, extended=T, recursive=T, force=T)
	sampleNames(rgset) <- pData(rgset)$Sample_Name
	message("Data read.")

	### QC raw data ###
	plotControlProbes(rgset, outdir=file.path(outdir, 'qc'))
	plotSimplifiedQC(preprocessRaw(rgset), prefix=file.path(outdir, 'qc', 'Raw'))
	
	### Genotype samples ###
	genotypeSamples(rgset, outdir=file.path(outdir, 'qc'))
	
	### Fetching array annotation ###
	array.type <- annotation(rgset)['array'] # IlluminaHumanMethylation450k / IlluminaHumanMethylationEPIC
	array.annot <- getAnnotation(rgset)
	array.annot.gr <- GRanges(array.annot$chr, ranges=IRanges(array.annot$pos-1, array.annot$pos), strand=NULL, name=array.annot$Name, score=0, array.annot[,! colnames(array.annot) %in% c('chr', 'pos', 'strand', 'Name')], seqinfo=Seqinfo(sort(unique(array.annot$chr))))
#	array.annot.gr <- array.annot.gr[,-match('Name', colnames(mcols(array.annot.gr)))]
#	array.annot.gr <- sort(array.annot.gr)
	annot.bed <- gr2bed(array.annot.gr)
	
	## Output annotation ##
	write.table(array.annot, file.path(outdir, 'annotation.txt'), sep="\t", quote=F, row.names=T)
	
	### Filter data ###
	
	## Flag low quality probes ##
	lowQ <- detectionP(rgset) > 0.01
	lowQ.probes <- rowSums(lowQ)/ncol(lowQ) >= 0.5
	lowQ.probes <- names(lowQ.probes)[lowQ.probes]
	annot.bed[lowQ.probes, 'score'] <- annot.bed[lowQ.probes, 'score'] + 8

	## Flag crossreactive probes ##
	crossreactive <- unlist(crossreactive[[array.type]])
	crossreactive <- crossreactive[crossreactive %in% annot.bed$name]
	annot.bed[crossreactive, 'score'] <- annot.bed[crossreactive, 'score'] + 2
	#o#
	
	## Flag blacklisted probes ##
	if (!is.null(blacklist)) {
		blacklisted.probes <- .parse_blacklist(blacklist)
		flag <- blacklisted.probes$flag
		flag <- flag[flag %in% rownames(annot.bed)]
		if (length(flag > 0)) {
			annot.bed[flag, 'score'] <- annot.bed[flag, 'score'] + 4
		}
		remove <- blacklisted.probes$remove
		remove <- remove[row.names(remove) %in% annot.bed$name,]
		remove <- remove[,colnames(remove) %in% targets$Sample_Name]
		if (nrow(remove) > 0) {
			remove <- remove[rownames(remove) %in% rownames(lowQ),]
			lowQ[row.names(remove), colnames(remove)] <- lowQ[row.names(remove), colnames(remove)] | remove
		}
		#o#
	}
	
	### Flag polymorphic probes ###
	if (removeSNPs) {
		af <- 1/ncol(rgset)
		polymorphic <- polymorphic[[array.type]]$CpG_SBE
		polymorphic <- !is.na(polymorphic[polymorphic[,paste(c(population, 'AF'), collapse="_")] > af,1])
		polymorphic <- polymorphic[polymorphic %in% annot.bed$name]
		annot.bed[polymorphic, 'score'] <- annot.bed[polymorphic, 'score'] + 1
	}
	
	### Add flag to array.annot.gr ###
	array.annot.gr$score <- annot.bed[names(array.annot.gr), 'score']

	### Sex determination ###
	ratio.meth <- mapToGenome(rgset, mergeManifest=T)
	gender <- getSex(ratio.meth)
	pData(rgset)$predictedSex <- factor(gender$predictedSex)

	### Reformat targets ###
	pdata <- pData(rgset)
	write.csv(pdata, file.path(outdir, 'extended_sample_sheet.csv'), quote=F, row.names=F)
	noisy.vars <- c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Pool_ID', 'Array', 'Basename', 'filenames')
	noisy.vars <- sapply(pdata, function(x) all(is.na(x)))
	noisy.vars <- c(colnames(pdata)[noisy.vars], 'Basename', 'filenames')
	if (!usePredictedSex) noisy.vars <- c(noisy.vars, 'predictedSex')
	targets <- pdata[,!colnames(pdata) %in% noisy.vars]
	targets$ArrayRow <- gsub('C..', '', pdata$Array)
	targets$ArrayColumn <- gsub('R..', '', pdata$Array)
	targets$Slide <- as.factor(targets$Slide)
	not.unique <- sapply(targets, function(x) length(unique(x)))
	not.unique <- not.unique > 1
	targets <- targets[, not.unique]
#	plot.vars <- unique(c('Slide', 'ArrayRow', 'ArrayColumn', batch.vars))

	### Output raw tables ###
	raw.betas <- getBeta(rgset)
	betas.sd <- apply(raw.betas[rownames(raw.betas) %in% rownames(annot.bed)[annot.bed$score==0],], 1, sd)
	betas.sd <- sort(betas.sd, decreasing=T, na.last=T)
	var.betas <- head(names(betas.sd), 50000)
	inspector(raw.betas, targets, subset=var.betas, outdir=file.path(outdir, 'qc'), prefix='Raw')

	annot.bed <- annot.bed[rownames(annot.bed) %in% rownames(raw.betas),]
	raw.betas.bed <- data.frame(annot.bed, raw.betas[rownames(annot.bed),], stringsAsFactors=F, check.names=F)	
	saveRDS(rgset, file=file.path(outdir, 'rgset.rds'))
	gz <- gzfile(file.path(outdir, 'raw_betas.bed.gz'), 'w', compression=9)
	write.table(raw.betas.bed, gz, sep="\t", quote=F, row.names=F)
	close(gz)

	#### Process Data ####
	
	### Remove background ###
	message("Correcting the data...")
	norm.meth <- preprocessNoob(rgset, dyeCorr=T, dyeMethod='single')
	norm.meth <- preprocessSWAN(rgset, norm.meth)
	saveRDS(norm.meth, file=file.path(outdir, 'normalized_meth.rds'))
	message("Data corrected.")

	### QC processed data ###
	plotSimplifiedQC(norm.meth, prefix=file.path(outdir, 'qc', 'Normalized'))

	### CNV analysis ###
	if (runCNV) {
		CNV(norm.meth, exclude=names(array.annot.gr)[array.annot.gr!=0], outdir=file.path(outdir, 'qc', 'CNV_report'))
	}
	
	### Mask probes ###
	message("Masking unreliable values...")
	filtered.norm.meth <- norm.meth
	
	## Remove measures with low detection ###
#	assayDataElement(filtered.norm.meth, 'Meth')[lowQ] <- NA
#	assayDataElement(filtered.norm.meth, 'Unmeth')[lowQ] <- NA
	message("Masking done.")
	
	### Output preprocessed data ###
	message("Writing output...")
	processed.betas <- getBeta(filtered.norm.meth)
	processed.betas[lowQ] <- NA
#	colnames(processed.betas) <- targets[colnames(processed.betas), 'Sample_Name']
	inspector(processed.betas, targets, subset=var.betas, outdir=file.path(outdir, 'qc'), prefix='Normalized')
	betas.bed <- data.frame(annot.bed, processed.betas[rownames(annot.bed),], stringsAsFactors=F, check.names=F)
	processed.mval <- getM(filtered.norm.meth)
	processed.mval[lowQ] <- NA
#	colnames(processed.mval) <- targets[colnames(processed.mval), 'Sample_Name']
	mval.bed <- data.frame(annot.bed, processed.mval[rownames(annot.bed),], stringsAsFactors=F, check.names=F)
	saveRDS(filtered.norm.meth, file=file.path(outdir, 'filtered_normalized_meth.rds'))
	gz <- gzfile(file.path(outdir, 'filtered_normalized_betas.bed.gz'), 'w', compression=9)
	write.table(betas.bed, gz, sep="\t", quote=F, row.names=F)
	close(gz)
	gz <- gzfile(file.path(outdir, 'filtered_normalized_M.bed.gz'), 'w', compression=9)
	write.table(mval.bed, gz, sep="\t", quote=F, row.names=F)
	close(gz)
	message("Data ready in the output folder.")

	remove.cols <- c('Slide', 'ArrayRow', 'ArrayColumn')
	remove.cols <- remove.cols[!remove.cols %in% batch.vars]
	targets <- targets[,-match(remove.cols, colnames(targets))]

	return(list(M=processed.mval, betas=processed.betas, targets=targets, array.annot.gr=array.annot.gr))
}

.parse_blacklist <- function(blacklist)
{
	# Load blacklist
	remove <- read.delim(blacklist, sep='\t', stringsAsFactors=F, colClasses=c('character', 'character'), header=F)
	colnames(remove) <- c('probeID', 'sample')
#	remove <- remove[remove$probeID %in% annot.bed$name,]
	flag <- unique(remove[remove$sample=='',1])
	remove <- remove[remove$sample!='',]
	if (nrow(remove) > 0) {
		remove$value <- TRUE
		library(reshape2)
		remove <- dcast(probeID ~ sample, data=remove, fill=FALSE, value.var='value')
		row.names(remove) <- remove$probeID
		remove <- remove[,-1]
		remove <- as.matrix(remove)
	}
	return(list(flag=flag, remove=remove))
	#o#
}
