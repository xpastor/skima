report_session <- function(outdir=getwd())
{
	source('http://bioconductor.org/biocLite.R')
	session.tex <- file.path(outdir, 'session.tex')
	write(paste0('\\documentclass{report}\n\\title{\'yapima\' sessionInfo}\n\n\\usepackage{hyperref}\n\n\\begin{document}\n\\section*{\\centerline{\'yapima\' sessionInfo}}\n\\center{\\today}\n\n\\begin{itemize}\\raggedright\n  \\item Bioconductor version ', biocVersion(), '\n\\end{itemize}'),file=session.tex)
	write(toLatex(sessionInfo()), file=session.tex, append=T)
	write('\n\\end{document}', file=session.tex, append=T)
	setwd(outdir)
	texi2pdf(session.tex, clean=T)
}
