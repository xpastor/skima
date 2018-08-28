library(devtools)

#wd <- '/home/pastor/projects/InfiniumHumanMethylationIssues'
#wd <- '/icgc/dkfzlsdf/analysis/B080/pastor/InfiniumHumanMethylationIssues'
wd <- getwd()

#devtools::create(wd)
#setwd(wd)
#dir.create('data')
#dir.create('data-raw')

#devtools::use_data()
#devtools::use_data_raw()
# 450k, Chen et al.
#IlluminaHumanMethylation450kPAPER <- list()
crossreactivePAPER <- list()
polymorphicPAPER <- list()
crossreactive <- list()
polymorphic <- list()

source('https://gist.githubusercontent.com/schaunwheeler/5825002/raw/a7e1844d2abcb6c51f7d2479d5d9f64a473cb50f/xlsxToR.r')
download.file('http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48639-non-specific-probes-Illumina450k.xlsx', destfile=file.path(wd, 'data-raw', '48639-non-specific-probes-Illumina450k.xlsx'))
crossreactivePAPER$IlluminaHumanMethylation450k <- xlsxToR(file.path(wd, 'data-raw', '48639-non-specific-probes-Illumina450k.xlsx'), keep_sheets=c('nonspecific cg probes', 'nonspecific ch probes'), header=T)
crossreactivePAPER$IlluminaHumanMethylation450k <- lapply(crossreactivePAPER$IlluminaHumanMethylation450k, function(x) x[,1:4])

download.file('http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48640-polymorphic-CpGs-Illumina450k.xlsx', destfile=file.path(wd, 'data-raw', '48640-polymorphic-CpGs-Illumina450k.xlsx'))
polymorphicPAPER$IlluminaHumanMethylation450k <- xlsxToR(file.path(wd, 'data-raw', '48640-polymorphic-CpGs-Illumina450k.xlsx'), keep_sheets=c('Polymorphic CpGs & SNPs at SBE', 'SNPs within probes'), header=T)

#save(IlluminaHumanMethylation450kPAPER, file=file.path(wd, 'data', 'IlluminaHumanMethylation450kPAPER.RData'))
#devtools::use_data(IlluminaHumanMethylation450kPAPER)

crossreactive$IlluminaHumanMethylation450k <- list(cg_probes=crossreactivePAPER$IlluminaHumanMethylation450k$`nonspecific cg`[,1], ch_probes=crossreactivePAPER$IlluminaHumanMethylation450k$`nonspecific ch`[,1])

polymorphic_cg <- polymorphicPAPER$IlluminaHumanMethylation450k$`Polymorphic CpGs & SNPs at SBE`
distance <- c('MAPINFO+0'=0, 'MAPINFO+1'=1, 'SBE'=NA)
polymorphic_cg$BASE_FROM_MAPINFO <- distance[polymorphic_cg$BASE_FROM_MAPINFO]
sbe <- is.na(polymorphic_cg$BASE_FROM_MAPINFO)
distance <- matrix(c(-1,2,0,1), nrow=2, dimnames=list(c('I','II'), c('F','R')), byrow=T)
typeI <- polymorphic_cg$INFINIUM == 'I'
R <- polymorphic_cg$STRAND == 'R'
polymorphic_cg$BASE_FROM_MAPINFO[sbe & typeI & R] <- 2
polymorphic_cg$BASE_FROM_MAPINFO[sbe & typeI & !R] <- -1
polymorphic_cg$BASE_FROM_MAPINFO[sbe & !typeI & R] <- 1
polymorphic_cg$BASE_FROM_MAPINFO[sbe & !typeI & !R] <- 0
polymorphic_cg$BASE_FROM_SBE <- polymorphic_cg$BASE_FROM_SBE == 0
colnames(polymorphic_cg)[c(6,7)] <- c('DIST_FROM_MAPINFO', 'SBE')
polymorphic_cg <- polymorphic_cg[,c('PROBE', 'CHR', 'MAPINFO', 'INFINIUM', 'STRAND', 'DIST_FROM_MAPINFO', 'SBE', 'SNP_CHR', 'SNP_POS', 'SNP_ID', 'REF', 'ALT', 'AF', 'AFR_AF', 'AMR_AF', 'ASN_AF', 'EUR_AF')]

polymorphic_probe <- polymorphicPAPER$IlluminaHumanMethylation450k$`SNPs within probes`
polymorphic_probe$BASE_FROM_MAPINFO <- as.integer(gsub('MAPINFO.', '', polymorphic_probe$BASE_FROM_MAPINFO))
colnames(polymorphic_probe)[6] <- 'DIST_FROM_MAPINFO'
polymorphic_probe <- polymorphic_probe[,c('PROBE', 'CHR', 'MAPINFO', 'INFINIUM', 'STRAND', 'DIST_FROM_MAPINFO', 'SNP_CHR', 'SNP_POS', 'SNP_ID', 'REF', 'ALT', 'AF', 'AFR_AF', 'AMR_AF', 'ASN_AF', 'EUR_AF')]

polymorphic$IlluminaHumanMethylation450k <- list(CpG_SBE=polymorphic_cg, probe=polymorphic_probe)

#IlluminaHumanMethylation450k <- list(crossreactive=crossreactive, polymorphic=polymorphic)
#save(IlluminaHumanMethylation450k, file=file.path(wd, 'data', 'IlluminaHumanMethylation450k.RData'))
#devtools::use_data(IlluminaHumanMethylation450k)

# EPIC, McCartney et al.

#IlluminaHumanMethylationEPICPAPER <- list()

download.file('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830/bin/mmc2.txt', destfile=file.path(wd, 'data-raw', 'mmc2.txt'))
download.file('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830/bin/mmc3.txt', destfile=file.path(wd, 'data-raw', 'mmc3.txt'))
download.file('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830/bin/mmc1.txt', destfile=file.path(wd, 'data-raw', 'mmc1.txt'))

crossreactivePAPER$IlluminaHumanMethylationEPIC <- list()
crossreactivePAPER$IlluminaHumanMethylationEPIC$mmc2 <- scan(file.path(wd, 'data-raw', 'mmc2.txt'), what='character')
crossreactivePAPER$IlluminaHumanMethylationEPIC$mmc3 <- scan(file.path(wd, 'data-raw', 'mmc3.txt'), what='character')
polymorphicPAPER$IlluminaHumanMethylationEPIC <- list()
polymorphicPAPER$IlluminaHumanMethylationEPIC$mmc1 <- read.delim(file.path(wd, 'data-raw', 'mmc1.txt'), sep="\t", stringsAsFactors=F, header=T)

#save(IlluminaHumanMethylationEPICPAPER, file=file.path(wd, 'data', 'IlluminaHumanMethylationEPICPAPER.RData'))
devtools::use_data(crossreactivePAPER, overwrite=T)
devtools::use_data(polymorphicPAPER, overwrite=T)

crossreactive$IlluminaHumanMethylationEPIC <- crossreactivePAPER$IlluminaHumanMethylationEPIC
names(crossreactive$IlluminaHumanMethylationEPIC) <- c('cg_probes', 'ch_probes')

polymorphic$IlluminaHumanMethylationEPIC <- list(CpG_SBE=polymorphicPAPER$IlluminaHumanMethylationEPIC$mmc1)
colnames(polymorphic$IlluminaHumanMethylationEPIC$CpG_SBE)[c(1,4,9)] <- c('PROBE', 'INFINIUM', 'SNP_POS')

#IlluminaHumanMethylationEPIC <- list(crossreactive=crossreactive, polymorphic=polymorphic)
#save(IlluminaHumanMethylationEPIC, file=file.path(wd, 'data', 'IlluminaHumanMethylationEPIC.RData'))
devtools::use_data(crossreactive, overwrite=T)
devtools::use_data(polymorphic, overwrite=T)

tools::add_datalist(wd, force=T)
