run3.prediction.train.test <-
function(genome.train.file, genome.test.file, models.file, max.regressors=1, chain.iterations=10000)
{
	if (missing(genome.train.file)) stop("Name of the genome train file must be provided")
	if (missing(genome.test.file)) stop("Name of the genome test file must be provided")
	if (missing(models.file)) stop("Name of the models file must be provided")
	
	input.name <- paste("tmp.input.file.tt", round(runif(1, min=0, max=10000)), ".txt", sep='')
	# -----------------------------------------------
	# Determine dimensions of the genome file.
	# source("get.data.dims.R") 
	d.train <- get.data.dims(genome.train.file)
	d.test <- get.data.dims(genome.test.file)

        if(d.train$ncols == 0 || d.test$ncols == 0)
        	return()

	# ------------------------------------------------
	# Create the input file for the Prediction TrainTest step.
	# 
	write(paste("NumberOfVariables = ", d.train$ncols, sep=""), file=input.name, append=FALSE, sep="")
	write(paste("TrainDataFile = ", genome.train.file, sep=""), file=input.name, append=TRUE, sep="")
	 write(paste("TrainSampleSize = ", d.train$nrows, sep=""), file=input.name, append=TRUE, sep="")

        write(paste("TestDataFile = ", genome.test.file, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("TestSampleSize = ", d.test$nrows, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("MaxRegressors = ", max.regressors, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ModelsFile = ", models.file, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ChainIterations = ", chain.iterations, sep=""), file=input.name, append=TRUE, sep="")

	# ----------------------------------------------
	# Call the C MOSS function for step 2.

	# dyn.load("traintest.so")
	
	file.pass <- c(input.name)
	try(.C("rtraintest", as.character(file.pass), PACKAGE="traintest"))
	
	try(system(paste("rm ", input.name, sep="")))	

        numm <- as.integer(d.train$ncols)
        numm <- numm - 1
        return(paste(genome.test.file, ".", numm, ".fitted", sep=""))
}

