pre6.discretize.batch <- function(dir.file, dir.out, prefix.train, prefix.test=prefix.train, key.train="", key.test="", ending.train=".txt", ending.test=ending.train, train.output.append="binary", test.output.append="testbinary", splits.min=0.01, splits.inc=0.01, splits.max=0.99) {
#
# For all the train and test files in the directory dir.file,
# this function finds the splits in diallelic SNPs that are
# represented as three category discrete variables. The SNPs are
# dichotomized based on presense of 0 vs. absence of 0.
#
# Example:
#
# dir.file: directory where files with genotype information can be found.
# dir.out: output directory to which resulting file should be saved.
# prefix.train: the beginning of the file name for the pedegree file (up until chrom number)
# prefix.test: beginning of the file name for .dat file
# key.train: any keyword in the name of the pedegree file that distinguishes it from others
# key.test: any keyword in the name of the .dat file that distinguishes it from others
# ending.train: MUST be the extension of the pedegree file, including the dot ".".
# ending.test: MUST be the extension of the .dat file, including the dot ".".
#
# Outputs:
#
# <file.ped>_num<ending.ped> - in dir.out directory, the resultant binary file:
#      the SNP columns + last columns (but no user IDs will be recorded).

if(missing(dir.file)) stop("Name of input directory with geno files must be provided.")
if(missing(dir.out)) stop("Name of output directory must be provided")
if(missing(prefix.train)) stop("Prefix of the genos file name must be provided.")

# TODO: remove this line:
#source("get.file.name.R")
#source("get.chrom.num.R")


# *******************************************
# 1. Obtain all train and test files
all.train <- get.file.name(dir=dir.file, prefix=prefix.train, key=key.train, ending=ending.train)
all.test <- get.file.name(dir=dir.file, prefix=prefix.test, key=key.test, ending=ending.test)

if(length(all.train) == 0 || length(all.test) == 0)
	return()


# *******************************************
# 3. Match ped and .dat and run the pre5.genos2numeric()
 
chroms.train <- get.chrom.num(all.train, prefix=prefix.train)
chroms.test <- get.chrom.num(all.test, prefix=prefix.test)

chroms.common <- sort(intersect(chroms.train, chroms.test))

i <- 1
while (i <= length(chroms.common)) {
	curr.chrom <- chroms.common[i]

	curr.train <- all.train[match(curr.chrom, chroms.train)]
	curr.test <- all.test[match(curr.chrom, chroms.test)]

	print(paste("Discretizing: ", curr.train, " + ", curr.test, sep=""))

	pre6.discretize(file.train=curr.train, file.test=curr.test, dir.file=dir.file, dir.out=dir.out, train.output.append=train.output.append, test.output.append=test.output.append, splits.min=splits.min, splits.inc=splits.inc, splits.max=splits.max)
        i <- i + 1
}

}






