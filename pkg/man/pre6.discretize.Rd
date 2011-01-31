\name{pre6.discretize}
\Rdversion{1.1}
\alias{pre6.discretize}
\title{
Discretize the SNP data
}
\description{
This function finds the splits in diallelic SNPs that are represented as three category discrete variables. The SNPs are dichotomized based on presense of 0 vs. absence of 0. 
}
\usage{
pre6.discretize(file.train, file.test=file.train, dir.file, dir.out=dir.file, 
train.output.append="binary", test.output.append="testbinary", splits.min=0.01, 
splits.inc=0.01, splits.max=0.99)
}
\arguments{
  \item{file.train}{
	The name of the train input file that contains the SNPs (represented as 1, 2, and 3) together with binary outcomes. The binary outcome is assumed to occupy the last column in the file.
}
  \item{file.test}{
        The name of the test input file that contains the SNPs (represented as 1, 2, and 3) together with binary outcomes. The binary outcome is assumed to occupy the last column in the file. Defaults to \code{\var{file.train}}, if the \code{\var{file.test}} name is not specified.
}
  \item{dir.file}{
	The name of directory/path name which contains the \code{\var{file.train}} and \code{\var{file.test}} files. 
}
  \item{dir.out}{
	The name of directory/path name to which output files should be saved. Defaults to same location as \code{\var{dir.file}}
}
  \item{train.output.append}{
	The suffix that should be appended to the end of \code{\var{file.train}} filename to create the filename of the output file. The resultant output filename will be of the format <dir.out>/<file.train>.<train.output.append>.txt, where <file.train> is the filename without extension.
}
  \item{test.output.append}{
        The suffix that should be appended to the end of \code{\var{file.test}} filename to create the filename of the output test file. The resultant output filename will be of similar format to that created with \code{\var{train.output.append}}.
}
  \item{splits.min}{
	The minimum parameters that define the quantiles that are used to identify the possible splits. Unnecessary parameter when geno values are encoded as 3 categories.
}
  \item{splits.inc}{
        The minimum parameters that define the quantiles that are used to identify the possible splits. Unnecessary parameter when geno values are encoded as 3 categories.
}
  \item{splits.max}{
        The minimum parameters that define the quantiles that are used to identify the possible splits. Unnecessary parameter when geno values are encoded as 3 categories.
}

}
\details{
	 This is the function for Step1 of the MOSS algorithm. It requires its input files to contain genotype values as 1, 2, and 3, and the last column should represent the disease status (1 for CASE and 0 for CONTROL). The diallelic SNPs can be represented as three category discrete variables, since a segregating SNP site has three possible genotypes: 0/0, 0/1, and 1/1, where 0 is the wild type and 1 is the mutant allele. This algorithm dichotomizes the SNPs as presence of 0 vs. absence of 0, or as presence of 1 vs absence of 1. Two output files will be produced for train and test data, and the names of the output files will be returned. This function will also copy all files in directory \code{\var{dir.file}} that end with ".dat" and ".fam", into output directory \code{\var{dir.out}}.


Also this function will write an output file ending with .scores, which contains 4 columns of score values that are used for deciding the split. The .scores file is for debugging purposes: column 1 is the score for encoding 1 as 0 and 2&3 as 1; column 2: encode 1&2 as 0 and 3 as 1; column 3: encode 1&3 as 0 and 2 as 1; column 4 is the highest score chosen. 
}
\value{
\item{out$train}{The name of the train output file.}
\item{out$test}{The name of the test output file.}

}
\author{
Laurent Briollais, Adrian Dobra, Olga Vesselova
}
\note{
This step requires LINUX's \command{wc} functionality.
}

\seealso{
\code{\link{pre5.genos2numeric}}, \code{\link{pre5.genos2numeric.batch}}, 
\code{\link{pre6.discretize.batch}}, \code{\link{pre7.merge.genos}}, 
\code{\link{run1.moss.regression}}
}
\examples{

# Split fname into 2 parts for easy deletion of files afterwards.
fstart <- "randvals" 
fname <- paste(fstart, ".txt", sep="")

# Create a random matrix of 1,2,3s, and last column of 0s and 1s:
n.rows <- 6
n.cols <- 8
m.part1 <- matrix(round(runif(n.rows*n.cols, min=1, max=3)), n.rows)
m.part2 <- round(runif(n.rows))
#write.table(cbind(m.part1, m.part2), file=fname, append=FALSE, sep="\t", 
#col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(m.part1, m.part2), file=fname, append=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

outname <- pre6.discretize(file.train=fname, dir.file=".")

try(system(paste("rm ", fstart, "*", sep="")))
}

\keyword{misc}
