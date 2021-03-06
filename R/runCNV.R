## Required variables: ##
## pipeline_dir
## qcdir
## exclude
## raw.betas
## filtered.norm.meth
## processed.betas
## pdata
## genotype.betas

#source(file.path(pipeline_dir, 'functions.R'))

#library(GenomicRanges)
#library(cluster)
#library(parallel)

### CNV analysis ###
CNV <- function(mset, exclude=NULL, outdir=getwd())
{
	message('Running CNV analysis...')
	library(conumee)
	library(minfi)
#	cnv.dir <- file.path(outdir, 'CNV_report')
	dir.create(outdir, recursive=T)
	cnv.intensity <- getMeth(mset) + getUnmeth(mset)
	colnames(cnv.intensity) <- paste(colnames(cnv.intensity), 'intensity', sep='.')
	array.type <- annotation(mset)['array'] # IlluminaHumanMethylation450k / IlluminaHumanMethylationEPIC
	if (array.type == 'IlluminaHumanMethylation450k') {
		library(CopyNumber450kData)
		data(RGcontrolSetEx)
		#annotation(filtered.norm.meth)$array == 'IlluminaHumanMethylationEPIC'
	} else if (array.type == 'IlluminaHumanMethylationEPIC') {
		library(GEOquery)
		geos <- c('GSM2309177', 'GSM2309178', 'GSM2309179', 'GSM2309172', 'GSM2309180', 'GSM2309181', 'GSM2309182', 'GSM2309183', 'GSM2309184')
		geo.idatDir <- file.path(outdir, 'idat')
		dir.create(geo.idatDir, recursive=T)
		sapply(geos, getGEOSuppFiles, makeDirectory=F, baseDir=geo.idatDir)
		#	sapply(geos, getGEOidat, makeDirectory=F, baseDir=geo.idatDir)
		geo.files <- list.files(geo.idatDir, full.names=T)
		sapply(geo.files, gunzip)
		geo.files <- list.files(geo.idatDir, pattern='_Red', full.names=T)
		RGcontrolSetEx <- read.metharray(geo.files, extended=T)
		unlink(geo.idatDir, recursive=T, force=T)
	}
#	controls.norm <- preprocessENmix(RGcontrolSetEx)
	controls.norm <- preprocessNoob(RGcontrolSetEx)
	controls.norm <- preprocessSWAN(RGcontrolSetEx, mSet=controls.norm)
	
	exclude <- unique(exclude)
	exclude.gr <- if (!is.null(exclude)) getLocations(mset, lociNames=exclude)
	
	cnv <- CNV.load(as.data.frame(cnv.intensity), names=sub('\\.intensity$', '', colnames(cnv.intensity)))
	cnv.controls <- CNV.load(controls.norm)
	conumee.array <- c(IlluminaHumanMethylation450k='450k', IlluminaHumanMethylationEPIC='EPIC')
	cnv.annot <- CNV.create_anno(exclude_regions=exclude.gr, chrXY=T, array_type=conumee.array[array.type])
	for (pid in names(cnv)) {
		cat(pid)
		cnv.analysis <- CNV.fit(cnv[pid], cnv.controls, cnv.annot)
		cnv.analysis <- CNV.bin(cnv.analysis)
		cnv.analysis <- CNV.detail(cnv.analysis)
		cnv.analysis <- CNV.segment(cnv.analysis)
		cnv.res <- CNV.write(cnv.analysis, what='segments')
		cnv.probes <- CNV.write(cnv.analysis, what='probes')
		#	dir.create(file.path(qcdir, 'CNV_report', pid))
		#	write.table(cnv.res, file.path(qcdir, 'CNV_report', pid, paste(pid, 'CNV_report.txt', sep='_')), sep="\t", quote=F, row.names=F)
		write.table(cnv.res, file.path(outdir, paste(pid, 'CNV_report.txt', sep='_')), sep="\t", quote=F, row.names=F)
		pdf(file.path(outdir, paste(pid, 'CNV_report.pdf', sep='_')), width=20, height=5)
		#	png(file.path(qcdir, 'CNV_report', pid, paste(pid, 'whole_genome.png', sep='_')), width=2000, height=500)
		CNV.genomeplot(cnv.analysis)
		#	dev.off()
		for (chr in row.names(cnv.analysis@anno@genome)) {
			#		png(file.path(qcdir, 'CNV_report', pid, paste(pid, '_', chr, '.png', sep='')), width=2000, height=500)
			CNV.genomeplot(cnv.analysis, chr=chr)
			#		dev.off()
		}
		dev.off()
	}
	message('Finished.')
}
