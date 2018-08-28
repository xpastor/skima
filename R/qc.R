#' Control probes plots.
#'
#' Produces dotplots of the control probes provided by Illumina.
#'
#' A series of dotplots that may show issues in different steps of the array production are produced.
#'
#' 'staining.png' reports the efficiency of the staining step in each channel. These controls are independent of the hybridization and extension step.
#'
#' 'extension.png' reports the extension efficiency in a hairpin probe, and are therefore sample independent.
#'
#' 'hybridization.png' reports the overall performance of the Infinium Assay using syntetic targets in three different concentrations instead of amplified DNA.  The three different intensities related to the three concentrations should be clearly identifiable.
#'
#' 'target_removal.png' reports the efficiency of the stripping step after the extension. All target removal controls result in low signal compared to the hybridization controls if targets are removed efficiently.
#'
#' 'bisulfite_conversionI.png' and 'bisulfite_conversionII.png' report the efficiency of bisulfite conversion of non-CpG cytosines in the genomic DNA. For type I probes, PM probes should give high signal, while MM probes should give a signal close to background.
#'
#' 'specificityI.png' and 'specificityII.png' report potential nonspecific primer extension and are designed against non polymorphic T sites. For type I probes, PM probes should give high signal, while MM probes should give a signal close to background.
#'
#' 'non_polymorphic.png' reports the overall performance of the assay by querying a particular base in a nonpolymorphic region of the genome, letting you compare assay performance across different samples.
#'
#' 'negative.png' reports the overall background of the signal. The detectin p-value for each probe is estimated based on the intensities of the negative control probes.
#'
#' See \url{https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf} for more details.
#'
#' @param rgset a RGChannelSet(Extended) object from 'minfi'
#' @param outdir path to output the plots
#'
#' @return None
#'
#' @export
#'
#' @importFrom minfi controlStripPlot
#' 
plotControlProbes <- function(rgset, outdir=getwd())
{
	library(minfi)
	dir.create(outdir, showWarnings=F, recursive=T)
	dev.width <- .get_dev_width(rgset)
	res <- 300
	png(file.path(outdir, 'staining.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='STAINING', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'extension.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='EXTENSION', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'hybridization.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='HYBRIDIZATION', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'target_removal.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='TARGET REMOVAL', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'bisulfite_conversionI.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='BISULFITE CONVERSION I', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'bisulfite_conversionII.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='BISULFITE CONVERSION II', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'specificityI.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='SPECIFICITY I', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'specificityII.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='SPECIFICITY II', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'non_polymorphic.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='NON-POLYMORPHIC', sampNames=pData(rgset)$Sample_Name)
	dev.off()
	png(file.path(outdir, 'negative.png'), height=dev.width, width=7, units='in', res=res)
	controlStripPlot(rgset, controls='NEGATIVE', sampNames=pData(rgset)$Sample_Name)
	dev.off()
}

#' Median intensities plot.
#'
#' Plots the chipwide medians of the Meth and Unmeth channels
#'
#' @param mset a MethylSet(Extended) object from package 'minfi'
#' @param prefix prefix of the filename to save the file
#'
#' @return None
#'
#' @export
plotSimplifiedQC <- function(mset, prefix)
{
	library(minfi)
	med.int <- DataFrame(mMed=apply(log2(getMeth(mset)), 2, median, na.rm=T), uMed=apply(log2(getUnmeth(mset)), 2, median, na.rm=T))
	png(paste0(c(prefix, 'median_intensities.png'), collapse='_'), width=7, height=7, units='in', res=300)
	plotQC(med.int)
	dev.off()
}

#' Per sample density plot
#'
#' @details
#'
#' In 'mat', the features are in the rows and the samples in the columns.
#'
#' If 'groups' is specified, it has to be of the same length and order as the columns in mat, and the density lines will be coloured according to it.
#'
#' If no filename is provided, the plot will be forwarded to the current device.
#'
#' @param mat numeric matrix with measures per feature and sample, with samples in the columns and features in the rows
#' @param groups vector with annotation for the samples
#' @param filename character. File to save the density plot
#' @param main title of the density plot
#'
#' @export
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
plotSamplesDistribution <- function(mat, groups=NULL, filename=NULL, main='Density')
{
	library(ggplot2)
	library(reshape2)
	if (!is.null(filename)) {
		outdir <- dirname(filename)
		dir.create(outdir, showWarnings=F, recursive=T)
	}
	if (!is.null(groups) & (ncol(mat) != length(groups))) {
		stop('\'groups\' must be of the same length and order as the columns in \'mat\'.')
	}
	df <- melt(mat, value.name='Beta')
	gg <- ggplot(df, aes(x=Beta, group=Var2)) + theme_bw() + theme(panel.grid=element_blank(), panel.border=element_blank(), plot.title=element_text(hjust=0.5), legend.title=element_blank()) + labs(x='Beta value', y=element_blank(), title=main)
	if (is.null(groups)) {
		gg <- gg + geom_density()
	} else {
		names(groups) <- colnames(mat)
		annot.cols <- .annotation_colors(data.frame(groups, stringsAsFactors=F, row.names=names(groups)))
		df$Group <- groups[df$Var2]
		gg <- gg + geom_density(aes(col=Group)) + scale_color_manual(values=annot.cols$groups)
	}
	if (is.null(filename)) {
		gg
	} else {
		ggsave(gg, file=filename)
	}
}

#' Density heatmap
#'
#' @details
#' 'mat' is a numeric matrix of beta values (or values in general in the range [0,1]), with features in rows and samples in columns.
#'
#' If no filename is provided, the plot will be forwarded to the current graphical device.
#'
#' @inheritParams plotSamplesDistribution
#'
#' @export
#'
#' @importFrom ComplexHeatmap densityHeatmap
#'
plotBetaDensityHeatmap <- function(mat, filename=NULL, main=NULL)
{
	message('Plotting density heatmaps of Beta values...')
	if (!is.null(filename)) {
		outdir <- dirname(filename)
		dir.create(outdir, showWarnings=F, recursive=T)
	}
	library(ComplexHeatmap)
	dev.width <- .get_dev_width(mat, name='Density')
	heatmap.params <- list(cluster_rows=F, col=c('lightyellow', 'black'), name='Density')
	if (!is.null(filename)) {
		png(filename, width=dev.width, units='in', height=7, res=300)
	}
	densityHeatmap(mat, ylab='Beta value', title=main, range=c(0,1), cluster_columns=T, show_column_dend=T, column_names_gp=gpar(fontsize=7))
	if (!is.null(filename)) {
		dev.off()
	}
	message('Finished.')
}

#' Plot samples correlation
#'
#' @details
#' Pearson correlation between all the samples is calculated.
#'
#' @param annot data frame with annotation to be plotted. 
#' @inheritParams plotSamplesDistribution
#'
#' @export
#'
plotSampleCorrelation <- function(mat, annot=NULL, filename=NULL, main=NULL)
{
	### Correlation between samples ###
	message('Plotting sample correlations...')
	if (!is.null(filename)) {
		outdir <- dirname(filename)
		dir.create(outdir, showWarnings=F, recursive=T)
	}
	sample.cor <- cor(mat, use='na.or.complete')
	dev.width <- .get_dev_width(sample.cor, name='Correlation', annotation_names=colnames(annot))
	dev.height <- .get_dev_width(sample.cor, name='AAA')
	if (!is.null(filename)) {
		png(filename, width=dev.width, height=dev.height, units='in', res=300)
	}
	annotatedHeatmap(mat=sample.cor, name='Correlation', column_annotation=annot, row_annotation=annot, column_title=main, col=c('lightyellow', 'black'))
	if (!is.null(filename)) {
		dev.off()
	}
	message('Finished.')
}

#' Data inspection plotting functions.
#'
#' @param subset a vector with feature identifiers to cluster samples on a limited number of features
#' @inheritParams pca_report
#'
#' @export
inspector <- function(mat, annot, subset=rownames(mat), outdir=getwd(), prefix=NULL)
{
#	library(ggplot2)
#	library(ComplexHeatmap)
	dir.create(outdir, showWarnings=F, recursive=T)
	plotSamplesDistribution(mat, filename=file.path(outdir, paste0(c(prefix, 'beta_distribution.png'), collapse='_')), main=prefix)
	plotBetaDensityHeatmap(mat, filename=file.path(outdir, paste0(c(prefix, 'beta_distribution_heatmap.png'), collapse='_')), main=prefix)
	plotSampleCorrelation(mat[subset,], annot=annot, filename=file.path(outdir, paste0(c(prefix, 'samples_correlation.png'), collapse='_')), main=prefix)
	pca_report(mat[subset,], annot=annot, outdir=outdir, prefix=prefix)
}

#' Sample genotyping.
#'
#' Reports the genotypes that can be inferred using the 'rs' probes.
#'
#' @inheritParams plotControlProbes
#'
#' @export
#'
#' @importFrom minfi getSnpBeta
#' @import biomaRt
#' @import fpc
#' @import cluster
#'
genotypeSamples <- function(rgset, outdir)
{
	message('Genotyping samples...')
	### Extract genotyping probes ###
	genotype.betas <- getSnpBeta(rgset)
	snps <- row.names(genotype.betas)
	library(biomaRt)
	#mart <- useMart('snp')
	#snp.db <- useMart('snp', dataset='hsapiens_snp')
	snp.db <- useMart('ENSEMBL_MART_SNP', dataset='hsapiens_snp', host="grch37.ensembl.org")
	snp.genotype <- getBM(c('chr_name', 'chrom_start', 'chrom_end', 'refsnp_id', 'allele'), filters='snp_filter', values=snps, mart=snp.db)
	rownames(snp.genotype) <- snp.genotype$refsnp_id
	snp.unmeth <- gsub('[GC/]', '', snp.genotype$allele)
	snp.meth <- ifelse(snp.unmeth=='A', 'G', 'C')
	ref <- gsub('/..*', '', snp.genotype$allele)
	genotypes.df <- data.frame(ref=ref, hipo=paste0(snp.unmeth, snp.unmeth), hemi=paste0(snp.meth,snp.unmeth), hyper=paste0(snp.meth,snp.meth), row.names=snp.genotype$refsnp_id, stringsAsFactors=F)
	genotype.sample <- function(snps.betas) {
	    snps <- row.names(snps.betas)
		library(fpc)
		library(cluster)
		asw <- numeric(3)
		for (k in 2:3) asw[[k]] <- pam(snps.betas,k) $ silinfo $ avg.width
		k.best <- which.max(asw)
		
		hc <- hclust(dist(snps.betas))
		haplotype <- cutree(hc, k.best)
		genotype.mean <- aggregate(snps.betas, list(haplotype), mean)
		hyper <- genotype.mean$Group.1[which.max(genotype.mean$x)]
		hipo <- genotype.mean$Group.1[which.min(genotype.mean$x)]
		haplotype <- ifelse(haplotype==hipo, 'hipo', ifelse(haplotype==hyper, 'hyper', 'hemi'))
	}
	smpl.genotypes <- apply(genotype.betas, 2, genotype.sample)
	snp.idx <- match(row.names(smpl.genotypes), row.names(genotypes.df))
	genotype.idx <- apply(smpl.genotypes, 2, match, colnames(genotypes.df))
	genotypes <- matrix(genotypes.df[snp.idx + (genotype.idx - 1)*nrow(genotypes.df)], nrow=nrow(smpl.genotypes), ncol=ncol(smpl.genotypes), dimnames=list(row.names(smpl.genotypes), colnames(smpl.genotypes)))
	genotypes <- data.frame(chrom=snp.genotype[rownames(genotypes), 'chr_name'], chromStart=snp.genotype[rownames(genotypes), 'chrom_start']-1, chromEnd=snp.genotype[rownames(genotypes), 'chrom_end'], name=row.names(genotypes), score=rep('.', nrow(genotypes)), strand=rep('.', nrow(genotypes)), genotypes, stringsAsFactors=F)
	colnames(genotypes)[1] <- paste0('#', colnames(genotypes)[1])
	write.table(genotypes, file.path(outdir, 'genotypes.bed'), sep="\t", row.names=F, quote=F)
	message('Finished.')
}
