write.citations <- function(array.type, diffMeth=F, outdir=getwd())
{
	biblib <- NULL

	### minfi package ###
	biblib <- citation('minfi')[1]

	### Probe filtering ###
	if (array.type == 'IlluminaHumanMethylation450k') {
		chen <- bibentry(bibtype='Article', author='Yi-an Chen, Mathieu Lemire, Sanaa Choufani, Darci T. Butcher, Daria Grafodatskayak, Brent W. Zanke, Steven Gallinger, Thomas J. Hudson and Rosanna Weksberg', title='Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium Human Methylation450 microarray.', journal='Epigenetics', year='2013', volume='8', number='2', pages='203-209')
		biblib <- c(biblib, chen)
	} else if (array.type == 'IlluminaHumanMethylationEPIC') {
		mccartney <- bibentry(bibtype='Article', author='Daniel L. McCartney, Rosie M. Walker, Stewart W. Morris, Andrew M. McIntosh, David J. Porteous and Kathryn L. Evans', title='Identification of polymorphic and off-target probe binding sites on the Illumina Infinium MethylationEPIC BeadChip', journal='Genomics Data', year='2016', volume='9', pages='22-24')
		biblib <- c(biblib, mccartney)
	}

	### ENmix background correction ###
#	biblib <- c(biblib, citation('ENmix'))

	### Noob background correction ###
	biblib <- c(biblib, citation('minfi')[4])

	### SWAN normalization ###
	biblib <- c(biblib, citation('minfi')[2])

	### EPiC arrays ###
	if (array.type == 'IlluminaHumanMethylationEPIC')
		biblib <- c(biblib, citation('minfi')[7])

	### Differential Methylation ###
	if (diffMeth) {
	## limma for standard differential methylation analysis ##
		biblib <- c(biblib, citation('limma'))

	## DMRcate for DMR detection
		biblib <- c(biblib, citation('DMRcate'))

	## Benjamini & Hochberg multiple testing correction (fdr) ##
		bh <- bibentry(bibtype='Article', author='Y Benjamini and Y Hochberg', title='Controlling the false discovery rate: a practical and powerful approach to multiple testing.', year='1995', journal='Journal of the Royal Statistical Society Series B', volume='57', pages='289-300')
		biblib <- c(biblib, bh)
	}

	### M-value/Beta-value conversion ###
	#if (diffMeth) {
		m2beta <- bibentry(bibtype='Article', author='Pan Du, Xiao Zhang, Chiang-Ching Huang, Nadereh Jafari, Warren A Kibbe, Lifang Hou and Simon M Lin', title='Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis.', year='2010', journal='BMC Bioinformatics', volume='11', pages='587')
		biblib <- c(biblib, m2beta)
	#}

	biblib.text <- format(biblib, 'text')
	writeLines(biblib.text, con=file.path(outdir, 'citations.txt'), sep="\n\n")
}
