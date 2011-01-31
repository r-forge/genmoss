\name{genos.clean.batch}
\Rdversion{1.1}
\alias{genos.clean.batch}
\title{
	Removes badly predicted SNPs by MaCH for all files
}
\description{
        For all files in \code{\var{dir.ped}}, does the same thing as \code{\link{pre5.genos2numeric.batch}}, only leaves genotypes the way they are, without categorizing them into 3 levels. Removes all SNPs that have missing or bad values. Intended to be done after imputation, to ensure consistency. Geno values should use letters A, T, C, G if \code{\var{letter.encoding}}=TRUE.
}
\usage{
genos.clean.batch(dir.ped, dir.dat = dir.ped, dir.out, prefix.ped, prefix.dat, 
key.ped = "", key.dat = "", ending.ped = ".txt", ending.dat = ".dat", 
num.nonsnp.col = 2, num.nonsnp.last.col = 1, letter.encoding = TRUE, 
save.ids.name = "patients.fam")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.ped}{
        The name of directory that contains all the .ped files.
}
  \item{dir.dat}{
        The name of directory that contains all the .dat files.
}
  \item{dir.out}{
        The name of output directory to which resulting files should be saved.
}
  \item{prefix.ped}{
        The beginning of the file name for the pedegree file (up until chrom number).
}
  \item{prefix.dat}{
        The beginning of the file name for .dat file (up until chrom number).
}
  \item{key.ped}{
        Any keyword in the name of the pedegree file that distinguishes it from other non-pedegree files.
}
  \item{key.dat}{
        Any keyword in the name of the .dat file that distinguishes it from others.
}
  \item{ending.ped}{
        MUST be the filename extension of the pedegree file, including the dot ".". For example, if your file is named "CGEM_2.txt", then set this variable to ".txt"; if your file is named "CGEM_2.ped", then set this variable to ".ped"; if your file is named "CGEM_2", then set this variable to "".
}
  \item{ending.dat}{
        MUST be the extension of the .dat file, including the dot ".".
}
  \item{num.nonsnp.col}{
        The number of leading columns that do not correspond to geno values. Ex. for MaCH1 input file format there are 5 non-snp columns; for MaCH1 output format .mlgeno it is 2; for Plink it is 6.
}
  \item{num.nonsnp.last.col}{
        The number of last columns that do not correspond to geno values. Ex. If last column is the disease status (0s and 1s), then set this variable to 1. If 2 last columns correspond to confounding variables, set the variable to 2.
}
  \item{letter.encoding}{
        Flag whether or not the ecoding used for Alleles is letters (A, C, T, G). If True, then does additional check for Alleles corresponding to the letters, and removes SNPs that contain any other symbols.
}
  \item{save.ids.name}{
	The name of the file to which all patient IDs should be saved.
}
}
\details{
	This function calls \code{\link{genos.clean}} for all the files in the directory, so that users do not have to call that function as many times as there are chromosomes.

	For all the .ped files that start with \code{\var{prefix.ped}}, contain \code{\var{key.ped}}, and end with \code{\var{ending.ped}} in the directory \code{\var{dir.ped}}; and for similarly obtained .dat files, this function removes all the SNPs that have not been properly imputed by MaCH, making sure that there are no missing/strange values. This function is needed since results of MaCH might contain weird symbols (like '2' can appear instead of A, T, C, G). This is only effective when \code{\var{letter.encoding}} = True. The reason for calling this function, and not \code{\link{pre5.genos2numeric}} is because you might wish to call other software packages on the fully imputed data, which will not need the data categorized into 3 levels.

Outputs the following files:
\preformatted{
 <file.ped>_clean<ending.ped> - in dir.out directory, the resultant file:
      the SNP columns + last columns (but no user IDs will be recorded).
 <file.dat>_clean.dat - in dir.out directory, the corresponding .dat file, will 
      be different from original <file.dat> if any bad SNPs get removed.
 <save.ids.name> - the patient IDs, if save.ids.name is not empty "".
}
}

\author{
Olia Vesselova
}

\seealso{
\code{\link{pre3.call.mach}}, \code{\link{pre5.genos2numeric}}, 
\code{\link{pre5.genos2numeric.batch}}
}
\examples{
print("See demo for pre5.genos2numeric()")
}
