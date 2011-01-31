run2.prediction.cvv <-
function(genome.file, models.file, max.regressors=1, cvv.fold=2, chain.iterations=10000)
{
	if (missing(genome.file)) stop("Name of the genome datafile must be provided")
	if (missing(models.file)) stop("Name of the  models file from MOSS (ending with .reg) must be provided")

	input.name <- paste("tmp.input.file.cvv.", round(runif(1, min=0, max=10000)), ".txt", sep='')
	# -----------------------------------------------
	# Determine dimensions of the genome file.
	# source("get.data.dims.R") 
	d <- get.data.dims(genome.file)
        if(d$ncols == 0)
        	return()	

	# ------------------------------------------------
	# Create the input file for the Prediction Cvv step.
	# 
	write(paste("NumberOfVariables = ", d$ncols, sep=""), file=input.name, append=FALSE, sep="")
	write(paste("DataFile = ", genome.file, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("SampleSize = ", d$nrows, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("MaxRegressors = ", max.regressors, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ModelsFile = ", models.file, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("CvvFold = ", cvv.fold, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ChainIterations = ", chain.iterations, sep=""), file=input.name, append=TRUE, sep="")

	# ----------------------------------------------
	# Call the C MOSS function for step 2.

	# dyn.load("cvv.so")
	
	file.pass <- c(input.name)
	try(.C("rcvv", as.character(file.pass), PACKAGE="cvv"))
	
	try(system(paste("rm ", input.name, sep="")))	

	numm <- as.integer(d$ncols)
	numm <- numm - 1
	return(paste(genome.file, ".", numm, ".cvv.txt", sep=""))

}

