#' Produces a report with annotated PCA plots.
#'
#' @details
#' For each variable in 'annot' a PCA plot is produced. All the PCA plots are the same, only the representation of the annotation changes.
#'
#' @inheritParams plotSampleCorrelation
#' @inheritParams plotSimplifiedQC
#' @inheritParams plotControlProbes
#' @param outdir path to output the plots; if 'NULL' then the PCA is printed in the current open device
#'
#' @export
#'
#' @import ggplot2
#' @importFrom reshape2 melt
pca_report <- function(mat, annot=NULL, outdir=NULL, prefix=NULL)
{
	### PCA analysis ###
	message('PCA plots...')
	if (!is.null(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
	mat.narm <- mat[! apply(is.na(mat), 1, any),]
	pca <- prcomp(t(mat.narm))
	n.comp <- min(ncol(pca$x), 6)
	if (n.comp > 2) {
		list.grobs <- lapply(seq(ncol(annot)), function(i) {
				groups <- annot[,i]
				names(groups) <- row.names(annot)
				annotatedPCA(pca, groups, colnames(annot)[i])
			})
		names(list.grobs) <- colnames(annot)
		if (!is.null(outdir)) {
			n.grobs <- ceiling(n.comp/2)
			width.unit <- 13/20
			width <- n.grobs*2*4 + width.unit
			for (plot.var in names(list.grobs))
				ggsave(filename=file.path(outdir, paste0(c(prefix, plot.var, 'PCA.png'), collapse='_')), plot=list.grobs[[plot.var]], width=width, scale=0.7)
		} else {
			print(marrangeGrob(list.grobs, ncol=2, nrow=3))
		}
	}
	if (!is.null(outdir)) {
		library(ggplot2)
		library(reshape2)
		categorical <- sapply(annot, class)
		categorical <- names(categorical)[categorical %in% c('character', 'factor')]
		importance <- summary(pca)$importance[2,]
		importance <- importance[importance>=0.01]
		pca.df <- as.data.frame(pca$x)
		pca.df <- pca.df[,names(importance), drop=F]
		pca.df <- cbind(pca.df, annot[rownames(pca.df),categorical,drop=F])
		#	pca.df$id <- rownames(pca.df)
		pca.df <- as.data.frame(pca.df)
		pca.df <- melt(pca.df, id.vars=categorical)
		pca.df <- melt(pca.df, id.vars=c('variable', 'value'))
		colnames(pca.df) <- c('PC', 'x', 'variable', 'value')
		pca.df$PC.text <- paste0(pca.df$PC, ', ', importance[pca.df$PC]*100, '%')
		ncomp <- 4
		pdf(file.path(outdir,(paste0(c(prefix, 'PCA.pdf'), collapse='_'))), width=10, height=7)
		for (i in seq(ceiling(length(importance)/ncomp))) {
			comp <- paste0('PC', (i*ncomp-(ncomp-1)):(i*ncomp))
			tmp <- pca.df[pca.df$PC %in% comp,]
			print(ggplot(tmp, aes(x=value, y=x)) + geom_boxplot() + facet_grid(PC.text ~ variable, scales='free_x') + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1), panel.grid=element_blank()))
		}
		dev.off()
	}
	message('Finished.')
#	ggplot(pca.df, aes(x=value, y=x)) + geom_boxplot() + facet_grid(PC ~ variable, scales='free_x')
}

#' Annotated PCA plots.
#'
#' Function to produce PCA plots annotated for one feature.
#'
#' 'annotatedPCA' produces a plot with the first (up to) 6 components, 2 for each scatterplot, coloured following the annotation provided in 'groups'. In the axis of each scatterplot the annotation provided is plotted as scatterplot or as density plot in relation to the coordinates of the annotated component.
#'
#' @param groups a vector with annotation for the features
#' @param main a string to be used as title for the plot
#'
#' @rdname annotatedPCA
#' @export
#'
#' @importFrom grid grid.newpage textGrob
#' @importFrom gridExtra arrangeGrob
setGeneric('annotatedPCA', function(...) {
	standardGeneric('annotatedPCA')
})

#' @param x a named list with slots 'x' and 'importance'. 'x' stores a matrix  with the value of the rotated data. 'importance' is a numeric vector with the contribution of each component with values in [0,1]. The number of elements in 'importance' must be equal to the number of columns in 'x'.
#'
#' @export
setMethod('annotatedPCA', 'list', function(x, groups, main=NULL, components=1:6, return.grob=FALSE, ...)
{
	pca <- x
	components <- components[components <= ncol(pca$x)]
	n.comp <- min(ncol(pca$x), length(components))
	if (n.comp > 1) {
		library(grid)
		library(gridExtra)
		grid.newpage()
		legend.plot <- .pca.grob(pca, groups, 1, 2, legend=T)
		main <- textGrob(main, gp=gpar(fontsize=20), just='top')
		n.grobs <- n.comp%/%2
		rep.comp <- n.comp < 6 & n.comp%%2 == 1
		list.grobs <- list()
		for (i in seq(n.grobs)) {
			list.grobs[[i]] <- .pca.grob(pca, groups, components[2*i-1], components[2*i])
		}
		if (rep.comp) {
			list.grobs[[length(list.grobs)+1]] <- .pca.grob(pca, groups, components[n.comp-1], components[n.comp])
		}
		widths <- rep(4, length(list.grobs))
		list.grobs[[length(list.grobs) + 1]] <- .legend.grob(legend.plot)
		widths <- c(widths, 1)
#		return(arrangeGrob(.pca.grob(pca, groups, 1, 2), .pca.grob(pca, groups, 3, 4), .pca.grob(pca, groups, 5, 6), .legend.grob(legend.plot), widths=c(4, 4, 4, 1), ncol=4, top=main))
		pca.grob <- arrangeGrob(grobs=list.grobs, widths=widths, ncol=length(widths), top=main)
		if (return.grob) {
			return(pca.grob)
		} else {
			plot(pca.grob)
		}
	}
})

#' @param obj an object of the specified class
#'
#' @export
setMethod('annotatedPCA', 'prcomp', function(obj, groups, main=NULL, components=1:6, return.grob=FALSE, ...)
{
	pca <- list(x = obj$x, importance = summary(obj)$importance[2,])
	annotatedPCA(pca, groups = groups, main = main, components = components, return.grob = return.grob, ...)
})


#' @export
setMethod('annotatedPCA', 'princomp', function(obj, groups, main=NULL, components=1:6, return.grob=FALSE, ...)
{
	importance <- obj$sdev^2/sum(obj$sdev^2)
	pca <- list(x = obj$scores, importance = importance)
	annotatedPCA(pca, groups = groups, main = main, components = components, return.grob = return.grob, ...)
})

#' @param mat a numeric matrix, with samples in columns and variables in rows.
#' @inheritParams stats::prcomp
#'
#' @export
setMethod('annotatedPCA', 'matrix', function(mat, groups, main=NULL, components=1:6, return.grob=FALSE, scale.=TRUE, center=TRUE, ...)
{
	mat.narm <- mat[! apply(is.na(mat), 1, any),]
	pca <- prcomp(t(mat.narm), scale.=scale., center=center)
	annotatedPCA(pca, groups = groups, main = main, components = components, return.grob = return.grob, ...)
})

#' @import ggplot2
#' @importFrom grid rectGrob
#' @importFrom gridExtra arrangeGrob
.pca.grob <- function(pca, groups, comp1=1, comp2=2, legend=F)
{
	library(ggplot2)
	library(gridExtra)
	xvals <- pca$x[,comp1]
	yvals <- pca$x[,comp2]
	importance <- pca$importance
	df <- data.frame(xvals=xvals, yvals=yvals, groups=groups)
	xtitle <- paste0('PC', comp1, ' (', round(importance[comp1]*100, digits=1), '% of variance)')
	ytitle <- paste0('PC', comp2, ' (', round(importance[comp2]*100, digits=1), '% of variance)')

	theme_pca <- theme_bw() + theme(legend.position='none', panel.grid=element_blank(), panel.border=element_rect(colour='black'))
	pca.plot <- ggplot(df) + geom_point(mapping=aes(x=xvals, y=yvals, colour=groups)) + theme_pca + theme(axis.text=element_blank(), axis.ticks=element_blank(), plot.margin=unit(c(0,0,1,1), 'lines')) + xlab(xtitle) + ylab(ytitle)
	if (legend) {
		return(legend.plot <- pca.plot + theme(legend.position='right', legend.key=element_blank()) + labs(colour=''))
	} else {
		theme_aux_plot <- theme_pca + theme(panel.border=element_blank(), axis.title.x=element_blank())
		xplot <- ggplot(df)
		yplot <- ggplot(df)
		if (class(groups) %in% c('factor', 'character')) {
			xplot <- xplot + geom_density(aes(x=xvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=xvals), na.rm=T, adjust=2, linetype='dotted') + theme_aux_plot + theme(axis.title.y=element_blank(), axis.text.y=element_text(colour='#00000000'), axis.text.x=element_blank(), axis.ticks=element_blank(), plot.margin=unit(c(2,-0.5,0,0), 'lines'))
			yplot <- yplot + geom_density(aes(x=yvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=yvals), na.rm=T, adjust=2, linetype='dotted') + theme_aux_plot + theme(axis.text.x=element_text(colour='#00000000'), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(-0.5,2,1,0), 'lines')) + coord_flip()
		} else {
			xplot <- xplot + geom_point(aes(x=xvals, y=groups), na.rm=T, size=0.5) + geom_smooth(aes(x=xvals, y=groups), na.rm=T, se=FALSE, method='loess', span=1) + theme_aux_plot + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=unit(c(2,0,0,0), 'lines')) + ylab('')
			yplot <- yplot + geom_point(aes(x=yvals, y=groups), na.rm=T, size=0.5) + geom_smooth(aes(x=yvals, y=groups), na.rm=T, se=FALSE, method='loess', span=1) + theme_aux_plot + theme(axis.title.y=element_blank(), axis.text.x=element_text(angle=270, vjust=0.5, hjust=0), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=unit(c(0,2,1,0), 'lines')) + coord_flip()
		}
		return(arrangeGrob(xplot, rectGrob(gp=gpar(lty='blank')), pca.plot, yplot, ncol=2, nrow=2, widths=c(3,1), heights=c(1,3), padding=unit(0, 'line')))
	}
}

#' @importFrom ggplot2 ggplot_gtable ggplot_build
.legend.grob <- function(a.gplot) {
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
	return(tmp$grobs[[leg]])
}
