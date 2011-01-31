\name{pre7.merge.genos}
\Rdversion{1.1}
\alias{pre7.merge.genos}
\title{
        Combine geno files across all chromosomes
}
\description{
Puts together all the genos files and their corresponding .dat files for all chromosomes. The files should have last column as the disease status, and the number of individuals (rows) must match across all files. Also the files are expected to have no leading non-snp columns. If they exist, they will be removed. The dat files are expected to have the SNP names in their second column. If the first column of .dat file is 'M', then it will be replaced by the chomosome number of the file name (the number that follows prefix.dat). This function tries to make sure that the geno files and dat files correspond.
}

\usage{
pre7.merge.genos(dir.file, dir.dat = dir.file, dir.out = dir.file, 
file.out = "CGEM_Breast_complete.txt", dat.out = "CGEM_Breast_complete.dat", 
prefix.file, prefix.dat, key.file = "", key.dat = "", ending.file = ".txt", 
ending.dat = ".dat", num.nonsnp.col = 0, num.nonsnp.last.col = 1, 
weak.check = FALSE, plan = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.file}{
        The name of directory containing files with geno information. The files in this directory must have their last column as the disease status.
}
  \item{dir.dat}{
        The name of directory containing .dat files. Should be a list of geno IDs, one ID per line, no header. Defaults to same directory as \code{\var{dir.genos}}.
}
  \item{dir.out}{
        The name of directory where the two output files will go. Defaults to same directory as \code{\var{dir.genos}}.
}
  \item{file.out}{
        The name of the output file which will contain the combined geno information and the last column will be the disease status.
}
  \item{dat.out}{
        The name of the output file which will contain all the corresponding SNP values.
}
  \item{prefix.file}{
        The string that appears at the beginning of all the geno input file names. The file names are expected to begin with \code{\var{prefix.file}}, and then be immediately followed by chromosome number, for example, in \code{\var{dir.file}} directory files named like :
\preformatted{
        "cgem_breast.21.pure.txt"
        "cgem_breast.5.pure.txt"
        "cgem_breast.24_and_25.txt"
         must have prefix="cgem_breast."
}
}
  \item{prefix.dat}{
        The string that appears at the beginning of all the .dat file names. Similarly to \code{\var{prefix.file}}, it must be immediately followed by the chromosome number.
}
  \item{key.file}{
        Any keyword in the name of the geno file that distinguishes it from other files.
}
  \item{key.dat}{
        Any keyword in the name of the .dat file that distinguishes it from other files.
}
  \item{ending.file}{
        The string with which all the geno filenames end.
}
  \item{ending.dat}{
        The string with which all the .dat filenames end.
}
  \item{num.nonsnp.col}{
        The number of leading columns in the .ped files that do not contain SNP values. The first columns of the file represent non-SNP values (like patient ID, gender, etc). For MaCH1 input format, the \code{\var{num.nonsnp.col=5}}, for PLINK it is 6 (due to extra disease status column).
}
  \item{num.nonsnp.last.col}{
        The number of last columns that do not correspond to geno values. Ex. If last column is the disease status (0s and 1s), then set this variable to 1. If 2 last columns correspond to confounding variables, set the variable to 2.
}
  \item{weak.check}{
        Since this function will try to check correspondence of the number of genos in the genos file to the .dat file, the function would expect there to be the same number of genos and .dat files. If you wish to by-pass these checks, set \code{\var{weak.check}}=TRUE, in which case only the total final number of the resultant geno and .dat files will be checked for consistency, and only a warning message will be printed if there is a problem.
}
  \item{plan}{
        Flag: if this option is TRUE, then this function will "do" nothing, but will simply print which files it plans to combine in which order, since combination step itself might take time for large files.
}
}
\details{
Puts together all the genos files and their corresponding .dat files for all chromosomes. The files should be tab separated and have last column as the disease status, and the number of individuals (rows) must match across all files. Also the files are expected to have no leading non-snp columns. If they exist, they will be removed. The dat files are expected to have the SNP names in their second column. If the first column of .dat file is 'M', then it will be replaced by the chomosome number of the file name (the number that follows prefix.dat). This function tries to make sure that the geno files and dat files correspond.


The resultant combined geno file will be saved into \code{\var{file.out}} and .dat file will be saved in \code{\var{dat.out}}.
}
\value{
	The FULL name of the combined result geno file (including the directory).
}
\note{
	The function makes use of LINUX commands: 'paste', 'cat', and 'wc'.
}
\author{
Olia Vesselova
}


\seealso{
\code{\link{pre6.discretize}}, \code{\link{pre6.discretize.batch}}, 
\code{\link{pre8.add.conf.var}}, \code{\link{pre9.split.train.test}}, 
\code{\link{pre9.split.train.test.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
