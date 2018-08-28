write.methods <- function(array.type='IlluminaHumanMethylation450k', blacklist='', af=NULL, population='mixed', diffMeth=F, outdir=getwd())
{	
	text <- "The Illumina Infinium methylation array analysis was done using the \'skima\' R package developed by Xavier Pastor (x.pastorhostench@dkfz-heidelberg.de) at DKFZ.\n\n"

### Preprocessing ###

	flag.ref <- NULL
	if (array.type == 'IlluminaHumanMethylation450k') {
		flag.ref <- '(Chen, 2013)'
	} else if (array.type == 'IlluminaHumanMethylationEPIC') {
		flag.ref <- '(McCartney, 2016)'
	}

	text <- paste0(text,
		'The array data were read into the R environment (R version ',
		R.version$major, '.', R.version$minor, ' \'', R.version$nickname, '\', ', R.version$year, '-', R.version$month, '-', R.version$day, ') using the \'minfi\' package ', .cite_package('minfi'), '.')

	text <- paste(text, 'The crossreactive probes', flag.ref, 'were flagged.')

	 if (!is.null(blacklist)) {
 		text <- paste0(text,
			'A list of blacklisted probes was also flagged (describe origin of blacklisted probes).')
	 }

	if (!is.null(af)) {
		text <- paste0(text,
#			' Probes with the single base extension (SBE) position overlaping SNPs with allele frequency higher than ', round(1/nrow(targets), digits=2),' in european populations ', flag.ref, ' were also flagged (n=', length(european.snps), ').')
			' Probes with the single base extension (SBE) position overlaping SNPs with allele frequency higher than ', round(af, digits=2),' in ', population,' populations ', flag.ref, '.')
	}

	text <- paste(text,
		'The intensities were adjusted with the \'noob\' background correction available in the \'minfi\' package with a default offset of 15, using a single sample approach and enabling dye bias correction (', .short_citation('minfi', 4), '). Normalization was done applying the \'SWAN\' normalization method (', .short_citation('minfi', 2), ').')
	#norm.meth <- preprocessENmix(rgset, bpParaEst='oob', dyeCorr=T, QCinfo=NULL, exQCsample=F, exQCcpg=F, exSample=NULL, exCpG=NULL, nCores=ncores)

	text <- paste(text,
		'Measures with a detection P-value higher than 0.01, as estimated by \'minfi\' and tipically of low quality, were masked. Probes with low quality in at least 50% of the samples were also flagged.')

	### Differential methylation analysis ###
	if (diffMeth) {
		limmaURL <- 'http://bioconductor.org/packages/3.1/bioc/vignettes/limma/inst/doc/usersguide.pdf'
		dmrcateURL <- 'https://www.bioconductor.org/packages/release/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf'
		text <- paste0(text, '\n',
			'The analysis of differentially methylated positions (DMP) was done on the M-values of the reliable probes (i.e. flag=0) using the biconductor \'limma\' package ', .cite_package('limma'), ' and multiple testing correction was applied (Benjamini, 1995) (protocol described in ', limmaURL, ').')
		text <- paste0(text, '\n',
			'The detection of differentially methylated regions (DMR) was done using the bioconductor \'DMRcate\' package ', .cite_package('DMRcate'), '. For two groups comparisons, the t statistics from the DMP analysis were used and the beta log fold change was computed running the standard \'limma\' workflow on the beta values. For comparisons with more than two groups the squared F statistics from the DMP analysis were provided and the beta log fold change was set to 0. All the other parameters were left as default.')
		text <- paste0(text, '\n\n',
			'When necessary, the beta-values were derived from the M-values as described by Du et. al. (Du, 2010).')
	}
#	return(outdir)
	write(text, file=file.path(outdir, 'methods.txt'))
}

.short_citation <- function(package, idx=1)
{
	cite <- citation(package)[idx]
	first <- cite$author[1]
	first <- gsub(' <.*', '', first)
	short_cite <- gsub('.* ', '', first)
	year <- cite$year
	if (! is.null(year)) short_cite <- paste(short_cite, year, sep=', ')
	return(short_cite)
}

.cite_package <- function(package, idx=1)
{
	short_cite <- .short_citation(package, idx)
	version <- paste0('v', packageVersion(package))
	cite_package <- paste0('(', short_cite, '; ', version, ')')
	return(cite_package)
}
