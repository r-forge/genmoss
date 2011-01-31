run1.moss.regression <-
function(genome.file, max.regressors=1, chain.iterations=10000, chain.replicates=5, cutoff.max=0.5, cutoff.min=0.001, prob.max=0.1, num.confounding.vars=0)
{
	if (missing(genome.file)) stop("Name of the genome datafile must be provided")

	input.name <- paste("tmp.input.file.", round(runif(1, min=0, max=10000)), ".txt", sep='')
	# -----------------------------------------------
	# Determine dimensions of the genome file.
	# source("get.data.dims.R") 
	d <- get.data.dims(genome.file)
	if(d$ncols == 0)
		return()

	# ------------------------------------------------
	# Hard code in the Shotgun parameters.
	shotgun.chain.replicates <- 1
	shotgun.cutoff.max <- 0.01
	shotgun.cutoff.min <- 0.0001
	shotgun.prob.max <- 0.1

	# ------------------------------------------------
	# Create the input file for MOSS algorithm.
	# 
	write(paste("DataFile = ", genome.file, sep=""), file=input.name, append=FALSE, sep="")	
	write(paste("NumberOfVariables = ", d$ncols ,sep=""), file=input.name, append=TRUE, sep="")
	write(paste("SampleSize = ", d$nrows, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("MaxRegressors = ", max.regressors, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ChainIterations = ", chain.iterations, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ChainReplicates = ", chain.replicates, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("CutoffMax = ", cutoff.max, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("CutoffMin = ", cutoff.min, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ProbMax = ", prob.max, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ShotgunChainReplicates = ", shotgun.chain.replicates, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ShotgunCutoffMax = ", shotgun.cutoff.max, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ShotgunCutoffMin = ", shotgun.cutoff.min, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ShotgunProbMax = ", shotgun.prob.max, sep=""), file=input.name, append=TRUE, sep="")
   write(paste("NumOfConfoundingVars = ", num.confounding.vars, sep=""), file=input.name, append=TRUE, sep="")


	# ----------------------------------------------
	# Call the C MOSS function for step 2.

	# dyn.load("shotgun.so")
	
	file.pass <- c(input.name)
	try(.C("rshotgun", as.character(file.pass), PACKAGE="shotgun"))
	
	try(system(paste("rm ", input.name, sep="")))	

	return(paste(genome.file, ".shotgun.", d$ncols, ".", max.regressors, ".reg", sep=""))
}

