.run_checks <- function(env=parent.frame())
{
	attach(env)
	if (is.null(idat_dir)) {
		stop('A directory with the IDAT files must be provided.')
	} else	if (!dir.exists(idat_dir)) {
		stop(paste0("\n\tThe directory with the IDAT files can not be accessed.\n\t", idat_dir))
	}
	
	if (is.null(sample.annotation)) {
		warning("\n\tNo sample sheet provided. The sample sheets present in the directory given at 'idat_dir' will be used.")
	} else if (!file.exists(sample.annotation)) {
		stop(sprintf("\n\tThe sample sheet does not exist.\n\t%s", sample.annotation))
	} else if (file.access(sample.annotation) == -1) {
		stop(sprintf("\n\tThe sample sheet could not be accessed.\n\t%s", sample.annotation))
	}
	
	if (!is.null(blacklist)) {
		if (! file.exists(blacklist) | file.access(blacklist) == -1) {
			stop(sprintf("\n\tThe file with blacklisted probes could not be accessed.\n\t%s", blacklist))
		}
	}
	
	if (!is.logical(runCNV)) {
		stop("\n\t'runCNV' must be a valid R boolean: T, F, TRUE or FALSE.")
	}
	
	if (!is.logical(diffMeth)) {
		stop("\n\t'diffMeth' must be a valid R boolean: T, F, TRUE or FALSE.")
	}
	
	if (!is.logical(removeSNPs)) {
		stop("\n\t'removeSNPs' must be a valid R boolean: T, F, TRUE or FALSE.")
	}
	if (removeSNPs) {
		ok <- sapply(polymorphic, function(x) grep('_AF$', colnames(x$CpG_SBE), value=T))
		ok <- sapply(ok, function(x) gsub('_AF$', '', x))
		all_ok <- unique(unlist(ok))
		exclusive <- sapply(ok, function(x) all_ok[!all_ok %in% x])
		exclusive <- melt(exclusive)
		rownames(exclusive) <- exclusive$value
		if (!(population %in% c('', all_ok))) {
			stop(paste0(population, ' is not a valid population. Specify one of ', paste0(c("''", all_ok), collapse=', '), '.'))
		} else if (population %in% rownames(exclusive)) {
			warning(paste0('Population \'', population, '\' not valid for ', exclusive[population, 'L1'], ' arrays.'))
		}
	}
	
	if (is.na(ncores)) {
		stop("\n\t'ncores' must be integer.")
	}
	
	if (is.na(seed)) {
		stop("\n\t'seed' must be integer.")
	}

	if (! dir.exists(outdir)) {
		success <- try(dir.create(outdir, recursive=T, mode='0770'), T)
		if (! success) {
			stop(paste0("\n\tThe output directory could not be created.\n\t", outdir))
		}
	} else {
		stop(paste0('\n\tThe directory \'', outdir, '\' already exists. Specify a non existing directory.'))
	}

	qcdir <- file.path(outdir, 'qc')
	dir.create(qcdir, recursive=T)

	detach(env)
}

.check_targets <- function(targets)
{
	categorical.vars <- colnames(targets)[sapply(targets, class) %in% c('character', 'factor')]
	no.replicates <- colnames(targets[,categorical.vars, drop=F])[apply(targets[,categorical.vars, drop=F], 2, function(x) any(table(x) == 1))]
	if (all(colnames(targets) %in% no.replicates)) {
		stop(paste0('The following variables contain groups without replicates: ', paste(no.replicates, collapse=', '), ". Replace them with 'NA' or reallocate them in other groups."))
	}
	return(no.replicates)
}
