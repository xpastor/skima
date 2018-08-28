#' Differential methylation analysis.
#'
#' Analysis of differentially methylated positions (DMPs) and detection of differentially methylated regions (DMRs).
#'
#' @param M matrix with M values
#' @param targets data frame with sample annotation to be used to run limma
#' @param array.annot.gr GRanges object with array annotation 
#' @param betas matrix with Beta values
#' @param interest.vars character vector with variables of interest, excluding the confounding factors
#' @param outdir path to save results
#' @param ncores integer. Number of cores to use to parallellize the DMR analysis
#' @param invisible a boolean to determine if the output of 'missMethyl::gometh' is returned
#' 
#' @inherit run_limma return
#'
#' @author Xavier Pastor \email{x.pastorhostench@dkfz.de}
#'
#' @export
#'
#' @importFrom minfi ilogit2
differentialMethylation <- function(M, targets, array.annot.gr, betas=ilogit2(M), interest.vars=NULL, outdir=getwd(), ncores=1, invisible=T)
{
	library(GenomicRanges)
	
	message('Starting differential methylation analysis...')
	
	dmp.dir <- file.path(outdir, 'DMP')
	dir.create(dmp.dir, recursive=T)
	dmr.dir <- file.path(outdir, 'DMR')
	dir.create(dmr.dir, recursive=T)
	
	# Limma analysis
	fit <- run_limma(M, targets)
	betafit <- run_limma(betas, targets)
	
	## Prepare data for DMR analysis
#	array.annot.gr <- sort(array.annot.gr[rownames(M)])
#	annot.bed <- gr2bed(array.annot.gr)

	if (is.null(interest.vars)) interest.vars <- colnames(targets)
	for (comparison in interest.vars) {
		message(comparison)
		var.coefs <- grep(paste0('^', comparison), colnames(fit$design), value=T)
		coef.name <- gsub(paste0('^', comparison), paste0(comparison, '.'), var.coefs)
		coef.name <- gsub('\\.$', '', coef.name)
		if (length(var.coefs) > 1) {
			coef.name <- comparison
		}
#		prefix <- file.path(dmp.dir, paste0('DMP_', coef.name))
		prefix <- file.path(dmp.dir, comparison)
		dir.create(prefix, recursive=T)
		prefix <- file.path(prefix, paste0('DMP_', coef.name))
		top <- report_significance(betas, fit, targets[,comparison], comparison, prefix, invisible=T)
		message('Writing significance table...')
		annot.bed <- gr2bed(array.annot.gr[rownames(top)])
		top <- cbind(annot.bed, top[rownames(annot.bed),], as.data.frame(mcols(array.annot.gr[rownames(annot.bed)])))
		top <- top[,!duplicated(colnames(top))]
		filename <- paste0(prefix, '.bed.gz')
		gz <- gzfile(filename, 'w', compression=9)
		write.table(top, gz, sep="\t", row.names=F, quote=F)
		close(gz)
		message('Finished')
		all.cpg <- top$adj.P.Val
		names(all.cpg) <- rownames(top)
		sig <- any(top$adj.P.Val <= 0.05)
		if (sig) {
			# Enrichment analysis
			report_enrichment(all.cpg, sig=0.05, collection='GO', prefix=paste0(prefix, '_'))
			report_enrichment(all.cpg, sig=0.05, collection='KEGG', prefix=paste0(prefix, '_'))
			# DMR analysis
			prefix <- file.path(dmr.dir, paste0('DMR_', coef.name))
			betatop <- topTable(betafit, coef=var.coefs, number=Inf, sort='none')
			dmr.gr <- DMR_detection(top, betatop, chr=as.character(seqnames(array.annot.gr[rownames(top)])), pos=end(array.annot.gr[rownames(top)]), ncores=ncores)
			report_DMRs(dmr.gr, array.annot.gr, betas=betas, annot=targets[,comparison,drop=F], prefix=prefix)
		}
	}
	message('Differential methylation analysis finished.')
	if (!invisible) return(fit)
}

#' Run limma on a dataset
#'
#' Generation of the fitted model with limma that can be use for any sort of differential analysis, as well as input for DMRcate
#'
#' The rows in 'targets' must be ordered like the columns in 'mat'. Before the fitted model is produced, 'targets' is checked for variables with absence of replicates and for colinearity.
#'
#' All the columns present in the 'targets' file will be included in the linear model.
#' 
#' @param mat numeric matrix
#' @param targets data frame with sample annotation
#'  
#' @return limma::eBayes
#'
#' @export
#'
#' @importFrom limma nonEstimable lmFit eBayes
#'
run_limma <- function(mat, targets, ...)
{
# Limma analysis
	library(limma)
	no.replicates <- .check_targets(targets)
	if (length(no.replicates) > 0) {
		warning(paste0('The following variables contain groups without replicates and will be removed before analysis: ', paste(no.replicates, collapse=', '), "."))
		targets <- targets[,! colnames(targets) %in% no.replicates]
	}
	my.formula <- paste(colnames(targets), collapse='+')
	my.formula <- as.formula(paste('~', my.formula))
	design <- model.matrix(my.formula, data=targets)
	n.e. <- nonEstimable(design)
	design <- design[, !colnames(design) %in% n.e.]
	mat <- mat[,rownames(design)]
	mat.narm <- mat[apply(!is.na(mat), 1, all),]
	fit <- lmFit(mat.narm, design, ...)
	fit <- eBayes(fit)
	return(fit)
}

#' Report significant results
#'
#' Significant results reported as PCA plots and dendrograms
#'
#' Extracts the significance analysis from a given variable. Features with an adjusted P-value lower or equal to 0.05 are taken as significant and used to summarize the results with PCA plots and hierarchical clustering.
#'
#' @param mat matrix used to report the results
#' @param fit fitLM object from limma
#' @param groups a named vector with the classification used to report results
#' @param coef a character with the tested condition to be reported
#' @param prefix character, the file prefix to report the results
#' @param invisible a boolean to determine if the output of 'missMethyl::gometh' is returned
#'
#' @return limma::topTable
#'
#' @export
#'
#' @importFrom limma topTable
#' @importFrom ggplot2 ggsave
#'
report_significance <- function(mat, fit, groups=NULL, coef=colnames(fit$design)[2], prefix=coef, invisible=T)
{
	library(ggplot2)
	var.coefs <- grep(paste0('^', coef), colnames(fit$design), value=T)
	top <- topTable(fit, coef=var.coefs, number=Inf, sort.by='none')
	coef.name <- gsub(paste0('^', coef), paste0(coef, '.'), var.coefs)
	if (length(var.coefs) == 1) {
		colnames(top) <- gsub('logFC', paste0(coef.name, '.logFC'), colnames(top))
	} else {
		colnames(top) <- gsub(paste0('^', coef), paste0(coef, '.'), colnames(top))
	}
	adjP <- grep('adj.P.Val', colnames(top))
	sig <- row.names(top)[top[,adjP] <= 0.05 & !is.na(top[,adjP])]
	if (length(sig) > 0) {
		sig.mat <- mat[sig,,drop=F]
		sig.mat <- sig.mat[apply(!is.na(sig.mat), 1, all),,drop=F]
		if (nrow(sig.mat) > 0) {
			pca <- prcomp(t(sig.mat), scale.=T, center=T)
#		groups <- targets[,coef]
#		names(groups) <- row.names(targets)
			ggsave(file=paste0(prefix,'_PCA.png'), annotatedPCA(pca, groups, coef), width=20)
			message('Reporting cluster...')
			scaled.mat <- apply(sig.mat, 1, scale, center=T)
			rownames(scaled.mat) <- colnames(mat)
			hc <- hclust(dist(scaled.mat))
			hc.width <- .get_dev_width(mat, name='', annotation_names=coef)
			png(paste0(prefix, '_cluster.png'), height=4, width=hc.width, units='in', res=300)
			hm_annot <- data.frame(groups, row.names=names(groups))
			colnames(hm_annot) <- coef
			annotatedHeatmap(matrix(nrow=0, ncol=ncol(mat), dimnames=list(NULL, colnames(mat))), column_annotation=hm_annot, cluster_columns=hc, column_dend_height=unit(5, 'cm'))
			dev.off()
			message('Finished.')
		}
	}
	return(top)
#	interest.vars <- .get_interest.vars(targets, c('Slide', batch.vars), usePredictedSex)
}

#' Report DMR analysis.
#'
#' Report the DMR tables in BED format and plots to inspect the results after DMR detection.
#'
#' DMRs are reported in BED format, sorted by chromosomal location.
#'
#' In order to produce compatible BED files, 'array.annot.gr' is sorted
#'
#' 'betas' is a matrix with the methylation level (Beta-value) of the CpGs. If 'betas' is provided, a heatmap is produced using the beta values corresponding to the CpG of the top DMRs, given the 'Stouffer', that sum up to 500 CpGs.
#'
#' If 'annot' is specified, the sample annotation is included in the heatmap.
#'
#' 'prefix' is the file prefix used to store the DMR results and plots.
#'
#' @param dmr.gr a GRanges object produced by 'extractRanges' from the 'DMRcate' package
#' @param array.annot.gr a GRanges object with array annotation
#' @param betas a matrix of beta values
#' @param annot a data frame with sample annotation
#' @param prefix character, file prefix to store the results
#'
#' @return None
#'
#' @export
#'
report_DMRs <- function(dmr.gr, array.annot.gr, betas=NULL, annot=NULL, prefix='DMR')
{
	library(ComplexHeatmap)
	dmr.bed <- gr2bed(dmr.gr)
	dmr.bed <- cbind(dmr.bed, as.data.frame(mcols(dmr.gr)))
	colnames(dmr.bed)[1] <- '#chrom'
	write.table(dmr.bed, paste0(prefix, '.bed'), sep='\t', quote=F, row.names=F)

	if (!is.null(betas)) {
		dmr.sig <- dmr.gr[order(dmr.gr$Stouffer)]
		L <- cumsum(dmr.sig$no.cpgs)
		n.dmr <- sum(L <= 500)
		dmr.sig <- head(dmr.sig, n.dmr)
		array.annot.gr <- array.annot.gr[rownames(betas),]
		array.annot.gr <- sort(array.annot.gr)
		overlap <- findOverlaps(array.annot.gr, dmr.sig)
		cpg.dmr <- array.annot.gr[queryHits(overlap),]
	
		hc.width <- .get_dev_width(betas, name='', annotation_names=colnames(annot))
		png(paste0(prefix, '.png'), height=10, width=hc.width, units='in', res=300)
		annotatedHeatmap(betas[names(cpg.dmr),], column_annotation=annot, row_annotation=data.frame(Chromosome=as.character(seqnames(cpg.dmr)), row.names=names(cpg.dmr), stringsAsFactors=F), cluster_rows=F, show_row_names=F, name='Beta', row_dend_side='left')
		dev.off()
	}
}

#' Report GO/KEGG analysis.
#'
#' @details
#' This function is a wrapper for the 'gometh' function from 'missMethyl' package that directly writes the results into a folder specified in 'prefix'.
#'
#' @seealso missMethyl::gometh
#'
#' @param all.cpg named vector, with significance or any value to subset significant probes
#' @param sig threshold used to subset the probes. Probes below or equal to the threshold are selected
#' @param collection 'GO' or 'KEGG'
#' @param prefix file prefix to write the results. It must point to a valid folder
#' @param invisible a boolean to determine if the output of 'missMethyl::gometh' is returned
#'
#' @inherit missMethyl::gometh return
#'
#' @export
#'
#' @importFrom missMethyl gometh
#'
report_enrichment <- function(all.cpg, sig=0.05, collection='GO', prefix=NULL, invisible=T)
{
	annot <- 'IlluminaHumanMethylation450kanno.ilmn12.hg19'
	library(annot, character.only=T)
	library(missMethyl)
	array.type <- NULL
	if (all(names(all.cpg) %in% rownames(Locations))) {
		array.type <- '450K'
	} else {
		array.type <- 'EPIC'
		annot <- 'IlluminaHumanMethylationEPICanno.ilm10b2.hg19'
		library(annot, character.only=T)
	}
	go <- gometh(names(all.cpg)[all.cpg<=sig], names(all.cpg), collection=collection, array.type=array.type, prior.prob=T)
	if (collection=='GO') {
#		load('/home/pastor/projects/packages/yapima/data/probe2GO.rda')
		terms <- unique(go$Ont)
		for (term in terms) {
			go.df <- go[go$Ont==term,]
			go.df <- go.df[order(go.df$FDR),]
#			GOdata <- new('topGOdata', ontology=term, allGenes=all.cpg, geneSel=function(x)x<=sig, nodeSize=1, annot=annFUN.gene2GO, gene2GO=probe2GO[[array.type]])
#			dag <- GOdata@graph
#			leaves.dag <- leaves(dag, degree.dir='in')
#			nodes.dag <- nodes(dag)
			go.df$GO_id <- rownames(go.df)
#			go.df$isLeaf <- F
#			go.df[rownames(go.df) %in% leaves.dag, 'isLeaf'] <- T
#			go.df[!rownames(go.df) %in% nodes.dag, 'isLeaf'] <- NA
			filename <- paste0(prefix, 'GO_', term, '.txt')
			write.table(go.df, filename, sep="\t", row.names=F, quote=F)
		}
	} else if (collection=='KEGG') {
		go$KEGG_id <- rownames(go)
		write.table(go, paste0(prefix, 'KEGG.txt'), sep="\t", quote=F, row.names=F)
	}
	if (!invisible) return(go)
}

reportGSA <- function()
{
	library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
	library(org.Hs.eg.db)
	library(limma)
	ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
	array.type <- if (all(all.cpg) %in% rownames(ann)) '450K' else 'EPIC'
	gsa <- gsa.meth(sig.cpg, all.cpg, collection=collection, array.type=array.type, prior.prob=T)
}
