#' Method DMR_detection.
#' @name DMR_detection
#' @rdname DMR_detection-methods
#' @exportMethod DMR_detection
setGeneric('DMR_detection', function(M, betas, ...){
	starndardGeneric('DMR_detection')
})

#' @rdname DMR_detection-methods
setMethod('DMR_detection', signature=c('data.frame', 'data.frame'),
	function(M, betas, chr, pos, significance=0.05, ncores=1) {
		library(DMRcate)
		R.utils::reassignInPackage('extractCoords', 'DMRcate', extractCoords)
		if (is.null(names(chr))) names(chr) <- rownames(M)
		if (is.null(names(pos))) names(pos) <- rownames(M)
		M <- M[!is.na(M$adj.P.Val),]
		betas <- betas[row.names(M),]
		chr <- chr[rownames(M)]
		pos <- pos[rownames(M)]
		betafc <- if('logFC' %in% colnames(betas)) {betas$logFC} else {0}
		stat <- if('t' %in% colnames(M)) {M$t} else {sqrt(M$F)}
		cpg.annot <- data.frame(ID=rownames(M), stat=stat, CHR=chr, pos=pos, betafc=betafc, indfdr=M$adj.P.Val, is.sig=M$adj.P.Val <= significance)
		class(cpg.annot) <- 'annot'
		dmr <- dmrcate(cpg.annot, mc.cores=ncores)
		dmr.gr <- extractRanges(dmr, genome='hg19')
		dmr.gr <- sort(dmr.gr)
		return(dmr.gr)
	}
)

setMethod('DMR_detection', signature=c('MArrayLM', 'MArrayLM'),
	function(M, betas, coef, chr, pos, significance=0.05, ncores=1) {
		library(limma)
		var.coefs <- grep(paste0('^', coef), colnames(M$design), value=T)
		topM <- topTable(fit, coef=var,coefs)
		topBeta <- topTable(betas, coef=var.coefs)
		dmr.gr <- DMR_detection(topM, topBeta, chr, pos, significance, ncores)
		return(dmr.gr)
	}
)
