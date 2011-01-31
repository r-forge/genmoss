pre6.discretize <- function(file.train, file.test=file.train, dir.file, dir.out=dir.file, train.output.append="binary", test.output.append="testbinary", splits.min=0.01, splits.inc=0.01, splits.max=0.99)
{
	if (missing(file.train)) stop("Name of the train datafile must be provided")
        if (missing(dir.file)) stop("The directory where train and test files are located must be provided (can even be an empty string: \"\", if file.train and file.test were given with full paths)")

        # -----------------------------------------------
	# Combine the full path input and output files
	full.train <- paste(dir.file, file.train, sep="/")
	full.test <- paste(dir.file, file.test, sep="/")

	# The output file name: remove the extension (last part after the last ".")
	splits.train <- unlist(strsplit(file.train, ".", fixed=T))
	start.train <- paste(splits.train[1:length(splits.train)-1], collapse=".")
	full.out.train <- paste(dir.out, "/", start.train, ".", train.output.append, ".txt", sep="")
	
	# The scores file name:
	full.out.scores <- paste(dir.out, "/", start.train, ".scores.txt", sep="")

        splits.test <- unlist(strsplit(file.test, ".", fixed=T))
        start.test <- paste(splits.test[1:length(splits.test)-1], collapse=".")
        full.out.test <- paste(dir.out, "/", start.test, ".", test.output.append, ".txt", sep="")
	
	input.name <- paste("tmp.step1.input.file.", round(runif(1, min=0, max=10000)), ".txt", sep='')
	# -----------------------------------------------
	# Determine dimensions of the train and test files.
	# source("get.data.dims.R") 
	d1 <- get.data.dims(full.train)
	if(d1$ncols == 0)
		return()
	d2 <- get.data.dims(full.test)
        if(d2$ncols == 0)
                return()

	# ------------------------------------------------
	# Create the input file for MOSS algorithm.
	# 
	write(paste("NumberOfVariables = ", d1$ncols ,sep=""), file=input.name, append=FALSE, sep="")
        write(paste("TrainDataFile = ", full.train ,sep=""), file=input.name, append=TRUE, sep="")
	write(paste("TrainSampleSize = ", d1$nrows, sep=""), file=input.name, append=TRUE, sep="")
        write(paste("TestDataFile = ", full.test ,sep=""), file=input.name, append=TRUE, sep="")
        write(paste("TestSampleSize = ", d2$nrows, sep=""), file=input.name, append=TRUE, sep="")

	write(paste("SplitsMin = ", splits.min, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("SplitsInc = ", splits.inc, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("SplitsMax = ", splits.max, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("TrainOutputFile = ", full.out.train, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("TestOutputFile = ", full.out.test, sep=""), file=input.name, append=TRUE, sep="")
	write(paste("ScoresFile = ", full.out.scores, sep=""), file=input.name, append=TRUE, sep="")



	# ----------------------------------------------
	# Call the C MakeSplits function for step 1.

	# dyn.load("ttsplits.so")
	
	file.pass <- c(input.name)
	try(.C("rsplits", as.character(file.pass), PACKAGE="ttsplits"))
	
	try(system(paste("rm ", input.name, sep="")))	

	# Copy over all .dat and .fam files
	get.file.copy(dir.in=dir.file, dir.out=dir.out, ending=".dat", verbal=FALSE)
	get.file.copy(dir.in=dir.file, dir.out=dir.out, ending=".fam", verbal=FALSE)

	return(list(train=full.out.train, test=full.out.test, score=full.out.scores))
}

