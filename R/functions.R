#' GenomicRanges to BED-like data frame.
#'
#' Function to convert GenomicRanges objects into BED-like data frames.
#'
#' The function produces a BED-like data frame where the 4th, 5th and 6th columns correspond to the name, the score and the strand associated to the feature. If the 'gr' object has a column called 'name', the name will be taken from this column, otherwise will be set to '.'. For the score, if a column named 'score' is present in the 'gr' object, the value will be taken from this column, otherwise it will be set to '.'. Any additional columns in the 'gr' object will be added as additional annotation from the 7th column on.
#'
#' @param gr a GenomicRanges object
#'
#' @return A data frame with at least six columns with the following info: chromosome, initial position, final position, name of the feature, score and strand.
#'
#' @export
#'
gr2bed <- function(gr)
{
	library(GenomicRanges)
	gr <- sort(gr)
	name <-	if ('name' %in% colnames(mcols(gr))) gr$name else '.'
	score <- if('score' %in% colnames(mcols(gr))) gr$score else '.'
	bed <- data.frame(chrom=as.character(seqnames(gr)), chromStart=start(gr)-1, chromEnd=end(gr), name=name, score=score, strand=strand(gr), stringsAsFactors=F)
	if ('name' %in% colnames(mcols(gr))) {
		rownames(bed) <- name
	}
	colnames(bed)[1] <- '#chrom'
	return(bed)
}

#' Guess array type from list of probes.
#'
#' Guesses the array type from a list of probes.
#'
#' @param x a character vector with a list of probes to determine the array
#'
#' @return A string with the array description as reported by 'minfi'.
#'
#' @export
#'
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#'
guessArrayType <- function(x)
{
	annot <- 'IlluminaHumanMethylation450kanno.ilmn12.hg19'
	library(annot, character.only=T)
	array.type <- ifelse(all(x %in% rownames(Locations)), 'IlluminaHumanMethylation450k', 'IlluminaHumanMethylationEPIC')
	return(array.type)
}
