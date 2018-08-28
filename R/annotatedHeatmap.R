#' Annotated heatmaps
#'
#' @param ... functions which define complex annotations or vectors of simple annotation. Values should be named arguments to be passed to 'Heatmap' function from 'ComplexHeatmap' package.
#' @param column_annotation,row_annotation a data frame. Each column will be treated as a simple annotation. The data frame must have column names.
#' @inheritParams ComplexHeatmap::Heatmap
#'
#' @export
#'
#' @import ComplexHeatmap
#' @import grid
annotatedHeatmap <- function(matrix, ..., column_annotation=NULL, row_annotation=NULL, column_names_gp=gpar(fontsize=7), row_names_gp=column_names_gp, row_dend_side='right', row_names_side='left', heatmap_legend_param=list(color_bar='continuous', legend_height=unit(3, 'cm')))
{
	library(ComplexHeatmap)
	heatmap.params <- list(matrix=matrix, column_names_gp=column_names_gp, row_names_gp=row_names_gp, row_dend_side=row_dend_side, row_names_side=row_names_side, heatmap_legend_param=heatmap_legend_param, ...)
	annotation.width <- unit(2, 'mm')
	if (! is.null(column_annotation)) {
		ha_cols <- .annotation_colors(column_annotation)
#		print(ha_cols)
		ha <- HeatmapAnnotation(df=column_annotation, col=ha_cols, gp=gpar(col='black'), na_col='white')
		heatmap.params$top_annotation <- ha
		annotation.width <- grobWidth(textGrob(colnames(column_annotation), gp=gpar(fontsize=12))) + unit(1, 'mm') - grobWidth(textGrob(rownames(matrix), gp=gpar(fontsize=7)))
	}
	hm <- do.call(Heatmap, args=heatmap.params)
	rows_hc <- ifelse(is.null(heatmap.params$cluster_rows), T, heatmap.params$cluster_rows)
	show_rows_hc <- ifelse(is.null(heatmap.params$show_row_dend), T, heatmap.params$show_row_dend)
	if (convertWidth(annotation.width, 'mm', valueOnly=T) < 0) annotation.width <- unit(2, 'mm')
	padding <- unit.c(unit(2, 'mm'), annotation.width, unit(2,'mm'), unit(2, 'mm'))
	if (is.null(row_annotation)) {
		draw(hm, padding=padding)
	} else {
		row_annot_params <- list(df=row_annotation, show_legend=!identical(row_annotation, column_annotation))
		if (identical(column_annotation, row_annotation)) {
#			col <- list()
			col <- ha_cols
#			for (anno in names(ha@anno_list)) {
#				if (class(row_annotation[,anno]) %in% c('character', 'factor')) {
#					col[[anno]] <- ha@anno_list[[anno]]@color_mapping@colors
#				} else {
#					col[[anno]] <- .annotation_colors(row_annotation[,anno, drop=F])[[anno]]
#				}
#			}
			row_annot_params <- c(row_annot_params, col=list(col), na_col='white')
#			row_annot_params$col <- col
			row_annot_params$gp <- gpar(lty='solid')
		} else {
			col <- .annotation_colors(row_annotation)
			row_annot_params <- c(row_annot_params, col=list(col))
		}
#		browser()
		row_annot <- do.call(rowAnnotation, row_annot_params)
		draw(hm + row_annot, row_dend_side=row_dend_side, padding=padding)
	}
	if (! is.null(column_annotation)) {
		for(ann in colnames(column_annotation)) {
			if(rows_hc & show_rows_hc) {
					decorate_annotation(ann, {grid.text(ann, unit(0, 'npc') - unit(2, 'mm'), just='right', gp=gpar(fontsize=10))})
			} else {
				  	decorate_annotation(ann, {grid.text(ann, unit(1, 'npc') + unit(2, 'mm'), just='left', gp=gpar(fontsize=10))})
			}
		}
	}
}

#' Annotated hierarchical clustering
#'
#' @param ... functions which define complex annotations or vectors of simple annotation. Values should be named arguments to be passed to 'Heatmap' function from 'ComplexHeatmap' package.
#' @param mat a matrix, with the samples organized in columns.
#' @param scale either a logical value or a numeric vector of length equal to the number of columns of ‘x’. Scale the data before clustering.
#' @param center either a logical value or a numeric vector of length equal to the number of columns of ‘x’. Center the data before clustering.
#' @param annotation a data frame. Each column will be treated as a simple annotation. The data frame must have column names.
#' @param dend_height height of the column cluster, should be a ‘unit’ object.
#' @inheritParams ComplexHeatmap::Heatmap
#'
#' @export
#'
#' @import ComplexHeatmap
#' @import grid
annotatedCluster <- function(mat, scale=T, center=T, annotation=NULL, dend_height = unit(5, 'cm'), silent=T, ...)
{
	scaled_mat <- apply(mat, 1, scale, scale=scale, center=center)
	rownames(scaled_mat) <- colnames(mat)
	hc <- hclust(dist(scaled_mat), ...)
	annotatedHeatmap(matrix(nrow = 0, ncol = nrow(scaled_mat), dimnames = list(NULL, rownames(scaled_mat))), column_annotation = annotation, cluster_columns = hc, column_dend_height = dend_height)
	if (!silent) return(hc)
}

#' @import grid
.get_dev_width <- function(mat, name='matrix_0', annotation_names=NULL, fontsize=7)
{
	library(grid)
	char.height <- convertHeight(grobHeight(textGrob("A", gp=gpar(fontsize = fontsize))), 'inch', valueOnly=T)
	width.mat <- ncol(mat) * char.height
	nchar.name <- max(nchar(unlist(strsplit(name, '\n'))))
	width.name <- convertWidth(grobHeight(textGrob('A', gp=gpar(fontsize=10))), 'inch', valueOnly=T) * nchar.name
	if (width.name < 0) width.name <- 0
	width.dev <- width.mat + convertWidth(unit(1, 'cm'), 'inch', valueOnly=T) + width.name + 4
	if (! is.null(annotation_names)) {
	    length.title <- max(nchar(annotation_names))
		width.title <- convertWidth(grobHeight(textGrob("A", gp=gpar(fontsize = 10))), 'inch', valueOnly=T) * length.title
		width.dev <- width.dev + width.title
	}
	return(width.dev)
}

.annotation_colors <- function(df)
{
#	library(circlize)
	cols <- c('darkgrey', 'black', 'red', 'yellow', 'blue', 'orange', 'cyan', 'magenta', 'darkgreen', 'khaki')
	annot_cols <- list()
	groups <- NULL
	categorical <- colnames(df)[sapply(df, class) %in% c('character', 'factor')]
	for (annot in categorical) {
		groups <- NULL
		x <- df[,annot]
		if (is.character(x)) {
			groups <- unique(x)
		} else if (is.factor(x)) {
			groups <- levels(x)
		}
		groups <- groups[!is.na(groups)]
		if (length(groups) <= length(cols)) {
			group_cols <- cols[seq(length(groups))]
		} else {
			group_cols <- rainbow(length(groups))
		}
		names(group_cols) <- groups
		annot_cols[[annot]] <- group_cols
	}
	other <- colnames(df)[!colnames(df) %in% categorical]
	cols <- cols[-1]
	numeric_cols <- seq_along(other)
	names(numeric_cols) <- other
#	numeric_cols <- lapply(numeric_cols, function(x) colorRampPalette(c('white', x)))
#	numeric_cols <- lapply(numeric_cols, function(x) colorRamp(c('white', x)))
#	numeric_cols <- lapply(numeric_cols, function(x) colorRamp2(c(0,10), c('white', x)))
	annot_cols <- c(annot_cols, numeric_cols)
	for (i in seq_along(other)) {
		my.col <- cols[(i-1) %% length(cols) + 1]
		my.col.func <- colorRampPalette(c('white', my.col))

		annot_cols[[other[i]]] <- .col_mapping(c('white', my.col), range(df[,other[i]], na.rm=T))
	}
	annot_cols <- annot_cols[colnames(df)]
	return(annot_cols)
}

.col_mapping <- function(colors, breaks)
{
	if (length(breaks) != 2) {
		stop('The vector \'ends\' must have 2 values.')
	}
	ends <- sort(breaks)
	attr <- list(breaks=breaks, colors=colors)
	fun <- function(x) {
		width.breaks <- breaks[2] - breaks[1]
		x[x < breaks[1]] <- breaks[1]
		x[x > breaks[2]] <- breaks[2]
		x <- (x-breaks[1])/width.breaks
		my.col <- colorRamp(colors, alpha=T)(x)
		my.col <- apply(my.col, 1, function(rgb) paste(as.hexmode(as.integer(rgb)), collapse=''))
		my.col <- paste0('#', toupper(my.col))
		my.col[is.na(x)] <- 'FFFFFFFF'
		return(my.col)
	}
	attributes(fun) <- attr
	return(fun)
}
