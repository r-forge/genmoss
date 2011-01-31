\name{pre5.genos2numeric}
\Rdversion{1.1}
\alias{pre5.genos2numeric}
\title{
        Categorize genotype data into 3 levels
}
\description{
	Categorizes genotype data into 3 levels, 1, 2, 3. Genos with two different Alleles are encoded as "2". Other genotypes are encoded as "1" or "3", where most frequent geno is "1". No missing values allowed, must be done after imputation. Geno values should use letters A, T, C, G if letter.encoding=TRUE. Also can work as a check for weird imputed values. 
}
\usage{
pre5.genos2numeric(file.ped, dir.ped, file.dat, dir.dat = dir.ped, dir.out, 
num.nonsnp.col = 2, num.nonsnp.last.col = 1, letter.encoding = TRUE, 
ped.has.ext = TRUE, dat.has.ext = TRUE, remove.bad.genos = FALSE, 
save.ids.name = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.ped}{
        The name of file with genotypes, after imputation.
}
  \item{dir.ped}{
        The name of directory where \code{\var{file.ped}} can be found.
}
  \item{file.dat}{
        The .dat file, should be tab separated, and no header.
}
  \item{dir.dat}{
        The name of directory where \code{\var{file.dat}} can be found. Defaults to \code{\var{dir.ped}}.
}
  \item{dir.out}{
        The name of output directory to which resulting file should be saved. The file will be named "Num.<file.ped>".
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
        
Categorizes genotype data into 3 levels, 1, 2, 3. Genos with two different Alleles are encoded as "2". Other genotypes are encoded as "1" or "3", where most frequent geno is "1". No missing values allowed, must be done after imputation. Geno values should use letters A, T, C, G if letter.encoding=TRUE. Also can work as a check for weird imputed values. For example, it is possible that an Allele is predicted by MaCH1 having value "2" (instead of A, T, C, or G) - it is best to remove SNPs that contain these weirdly imputed values.

The following files will be produced:
\preformatted{
 - <file.ped>_num<ending.ped> - in \code{\var{dir.out}} directory, the resultant 
    binary file: the SNP columns + last columns (but no user IDs will be 
    recorded), where <ending.ped> is the filename extension of file.ped.
 - <file.dat>_num.dat - in dir.out directory, the corresponding .dat file, will 
    be different from original <file.dat> if remove.bad.genos=TRUE.
 - <save.ids.name> - the patient IDs, if save.ids.name is not empty "".
}


}
\value{
        <file.ped>_num<ending.ped> filename - the name of the output file.
}
\author{
Olia Vesselova
}
\note{
 Note: in case of any bad values in the \code{\var{file.ped}} (ex. "NA", "0/0", "0", "1 1", etc), the output file Num_<file.ped> will still be produced, with '2' encoded by default in the place of bad input values, if \code{\var{remove.bad.genos}}=FALSE. Warning messages will be printed. If \code{\var{remove.bad.genos}}=TRUE, then these SNPs will be entirely removed, along with their names in the .dat file.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pre4.combine.case.control}}, 
\code{\link{pre4.combine.case.control.batch}}, 
\code{\link{pre5.genos2numeric.batch}}, 
\code{\link{pre6.discretize}}, \code{\link{pre6.discretize.batch}}
}

\examples{
print("See the demo 'gendemo'.")
}
