run4.save.prediction <- function(filename, outfile)
{
	# Returns the list with the following parameters:
	# tpr, fpr, rocarea
	# default value for output file is <filename>.plot.pdf

	if (missing(outfile)) outfile <- paste(filename,".plot.pdf",sep="")

        pdf(file=outfile) # prepare for saving
	
	outt <- run4.show.prediction(filename)

	dev.off()

	return(outt)
}

