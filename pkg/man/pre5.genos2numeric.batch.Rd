\name{pre5.genos2numeric.batch}
\Rdversion{1.1}
\alias{pre5.genos2numeric.batch}
\title{
        Categorize genotype data into 3 levels for each file
}
\description{
        For each .ped file in \code{\var{dir.ped}}, categorizes genotype data into 3 levels, 1, 2, 3. Genos with two different Alleles are encoded as "2". Other genotypes are encoded as "1" or "3", where most frequent geno is "1". No missing values allowed, must be done after imputation. Geno values should use letters A, T, C, G if letter.encoding=TRUE.
}
\usage{
pre5.genos2numeric.batch(dir.ped, dir.dat = dir.ped, dir.out, prefix.ped, 
prefix.dat, key.ped = "", key.dat = "", ending.ped = ".txt", ending.dat = ".dat", 
num.nonsnp.col = 2, num.nonsnp.last.col = 1, letter.encoding = TRUE, 
ped.has.ext = TRUE, dat.has.ext = TRUE, remove.bad.genos = FALSE, 
save.ids.name = "patients.fam")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.ped}{
        The name of directory where .ped files can be found.
}
  \item{dir.dat}{
        The name of directory where .dat files can be found.
}
  \item{dir.out}{
        The name of directory to which output files should go.
}
  \item{prefix.ped}{
        The beginning of the file name for the .ped pedegree file (up until chrom number).
}
  \item{prefix.dat}{
        The beginning of the file name for the .dat file (up until chrom number).
}
  \item{key.ped}{
        Any keyword in the name of the pedegree file that distinguishes it from other files.
}
  \item{key.dat}{
        Any keyword in the name of the .dat file that distinguishes it from other files.
}
  \item{ending.ped}{
        The ending of the pedegree filename.
}
  \item{ending.dat}{
        The ending of the .dat filename.
}
  \item{num.nonsnp.col}{
        The number of leading columns in the .ped files that do not contain SNP values. The first columns of the file represent non-SNP values (like patient ID, gender, etc). For MaCH1 input format, the \code{\var{num.nonsnp.col=5}}, for PLINK it is 6 (due to extra disease status column).
}
  \item{num.nonsnp.last.col}{
        The number of last columns that do not correspond to geno values. Ex. If last column is the disease status (0s and 1s), then set this variable to 1. If 2 last columns correspond to confounding variables, set the variable to 2.
}
  \item{letter.encoding}{
        Flag whether or not the ecoding used for Alleles is letters (A, C, T, G). If True, then does additional check for Alleles corresponding to the letters, and prints out warning messages if other symbols appear instead.
}
  \item{ped.has.ext}{
        Flag whether or not \code{\var{file.ped}} name has a filename extension (ex. ".ped", ".txt"). This is necessary for naming the output file.
}
  \item{dat.has.ext}{
        Flag whether or not file.dat name has a filename extension (ex. ".dat", ".txt").
}
  \item{remove.bad.genos}{
        Flag whether or not you want to remove a geno if at least one of its values is not valid (ex. "2" when only letters are expected, or "NA", etc). Warning: set this to TRUE only if the CASE and CONTROLs have been merged into the \code{\var{file.ped}}, (otherwise we do not want to remove some SNPs from CASE but not from CONTROL and generate two different .dat files).
}
  \item{save.ids.name}{
        The file name to which patient IDs should be saved. If not empty, then will save IDs of patients into another file with this name. Since dataset is generally split across many files, one chromosome each, the patient IDs should be the same across these files, thus it is enough to extract the patient ID ONCE, when running this code on the smallest chromosome. For runs on all other chromosomes, leave save.ids.name="" to save time and avoid redundant work. Could name output file as "patients.fam".
}
}
\details{
        For every pair of .dat and .ped files, categorizes genotype data into 3 levels, 1, 2, 3. Genos with two different Alleles are encoded as "2". Other genotypes are encoded as "1" or "3", where most frequent geno is "1". No missing values allowed, must be done after imputation. Geno values should use letters A, T, C, G if letter.encoding=TRUE. Also can work as a check for weird imputed values. For example, it is possible that an Allele is predicted by MaCH1 having value "2" (instead of A, T, C, or G) - it is best to remove SNPs that contain these weirdly imputed values.

The following files will be produced for each chromosome in the directory \code{\var{dir.ped}}:
\preformatted{
 - <file.ped>_num<ending.ped> - in \code{\var{dir.out}} directory, the resultant 
    binary file: the SNP columns + last columns (but no user IDs will be 
    recorded), where <ending.ped> is the filename extension of file.ped.
 - <file.dat>_num.dat - in dir.out directory, the corresponding .dat file, 
    will be different from original <file.dat> if remove.bad.genos=TRUE.
 - <save.ids.name> - the patient IDs, if save.ids.name is not empty "".
}

}
\author{
Olia Vesselova
}

\seealso{
\code{\link{pre4.combine.case.control}}, 
\code{\link{pre4.combine.case.control.batch}}, 
\code{\link{pre5.genos2numeric}}, \code{\link{pre6.discretize}}, 
\code{\link{pre6.discretize.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
