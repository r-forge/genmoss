\name{pre6.discretize.batch}
\Rdversion{1.1}
\alias{pre6.discretize.batch}
\title{
        Discretize the SNP data for all files
}
\description{
        For all the files, this function finds the splits in diallelic SNPs that are represented as three category discrete variables. The SNPs are dichotomized based on presense of 0 vs. absence of 0.
}
\usage{
pre6.discretize.batch(dir.file, dir.out, prefix.train, prefix.test = prefix.train, 
key.train = "", key.test = "", ending.train = ".txt", ending.test = ending.train, 
train.output.append = "binary", test.output.append = "testbinary", 
splits.min = 0.01, splits.inc = 0.01, splits.max = 0.99)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.file}{
        The name of directory that contains the train and test genotype files.
}
  \item{dir.out}{
        The name of directory to which output files should go.
}
  \item{prefix.train}{
        The beginning of the file name for the TRAIN file (up until chrom number).
}
  \item{prefix.test}{
        The beginning of the file name for the TEST file (up until chrom number).
}
  \item{key.train}{
        Any keyword in the name of the TRAIN file that distinguishes it from other files.
}
  \item{key.test}{
        Any keyword in the name of the TEST file that distinguishes it from other files.
}
  \item{ending.train}{
        MUST be the filename extension of the TRAIN file, including the dot ".". For example, if your file is named "CGEM_2.txt", then set this variable to ".txt"; if your file is named "CGEM_2.ped", then set this variable to ".ped".
}
  \item{ending.test}{
        MUST be the filename extension of the TEST file, including the dot ".". For example, if your file is named "CGEM_2.txt", then set this variable to ".txt"; if your file is named "CGEM_2.ped", then set this variable to ".ped".
}
  \item{train.output.append}{
        The suffix that should be appended to the end of TRAIN filename to create the filename of the output file. The resultant output filename will be of the format <dir.out>/<file.train>.<train.output.append>.txt, where <file.train> is the filename without extension.
}
  \item{test.output.append}{
        The suffix that should be appended to the end of TEST filename to create the filename of the output test file. The resultant output filename will be of similar format to that created with \code{\var{train.output.append}}.
}
  \item{splits.min}{
        The minimum parameters that define the quantiles that are used to identify the possible splits.
}
  \item{splits.inc}{
        The minimum parameters that define the quantiles that are used to identify the possible splits.
}
  \item{splits.max}{
        The minimum parameters that define the quantiles that are used to identify the possible splits.
}
}
\details{
        For all the TRAIN and TEST files, this is the function for Step1 of the MOSS algorithm. It requires its input files to contain genotype values as 1, 2, and 3, and the last column should represent the disease status (1 for CASE and 0 for CONTROL). The diallelic SNPs can be represented as three category discrete variables, since a segregating SNP site has three possible genotypes: 0/0, 0/1, and 1/1, where 0 is the wild type and 1 is the mutant allele. This algorithm dichotomizes the SNPs as presence of 0 vs. absence of 0, or as presence of 1 vs absence of 1. Two output files will be produced for train and test data.

This function will also copy all files in directory \code{\var{dir.file}} that end with ".dat" and ".fam", into output directory \code{\var{dir.out}}.

Since TRAIN and TEST set might not yet be split, it is all right to use the same file for both TRAIN and TEST.

}
\author{
Laurent Briollais, Adrian Dobra, Olga Vesselova
}


\seealso{
\code{\link{pre5.genos2numeric}}, \code{\link{pre5.genos2numeric.batch}}, 
\code{\link{pre6.discretize}}, \code{\link{pre7.merge.genos}}, 
\code{\link{run1.moss.regression}}

}
\examples{
print("See the demo 'gendemo'.")
}
