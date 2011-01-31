\name{genos.clean}
\Rdversion{1.1}
\alias{genos.clean}
\title{
	Removes badly predicted SNPs by MaCH
}
\description{
Same thing as \code{\link{pre5.genos2numeric}}, only leaves genotypes the way they are, without categorizing them into 3 levels. Removes all SNPs that have missing or bad values. Intended to be done after imputation, to ensure consistency. Geno values should use letters A, T, C, G if \code{\var{letter.encoding}}=TRUE.

}
\usage{
genos.clean(file.ped, ending.ped = ".txt", dir.ped, file.dat, ending.dat = ".dat", 
dir.dat = dir.ped, dir.out, num.nonsnp.col = 2, num.nonsnp.last.col = 1, 
letter.encoding = TRUE, save.ids.name = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.ped}{
        The name of file with genotypes, after imputation. Entries should be either tab or space separated.
}
  \item{ending.ped}{
        The extension of the \code{\var{file.ped}}, should contain the dot '.', if file has no ending, use an empty string "". This is needed to name the output file as <file.ped>_num<ending.ped>, where \code{\var{file.ped}} is without ending.
}
  \item{dir.ped}{
        The name of the directory where \code{\var{file.ped}} can be found.
}
  \item{file.dat}{
        The name of .dat file. This file should be tab separated, and no header.
}
  \item{ending.dat}{
        The extension of the \code{\var{file.dat}}, should contain the dot '.'. This is needed to name the output file as <file.dat>_num<ending.dat>, where \code{\var{file.dat}} is without ending.
}
  \item{dir.dat}{
        The name of directory where \code{\var{file.dat}} can be found. Defaults to \code{\var{dir.ped}}.
}
  \item{dir.out}{
        The name of output directory to which resulting files should be saved.
}
  \item{num.nonsnp.col}{
        The number of leading columns that do not correspond to geno values. Ex. for MaCH1 input file format there are 5 non-snp columns; for MaCH1 output format .mlgeno it is 2; for Plink it is 6.
}
  \item{num.nonsnp.last.col}{
        The number of last columns that do not correspond to geno values. Ex. If last column is the disease status (0s and 1s), then set this variable to 1. If 2 last columns correspond to confounding variables, set the variable to 2.
}
  \item{letter.encoding}{
        Flag whether or not the encoding used for Alleles is letters (A, C, T, G). If True, then does additional check for Alleles corresponding to the letters, and removes SNPs that contain values other than these 4 letters. Useful to eliminate 2s that may appear after MaCH1 imputation. 
}
  \item{save.ids.name}{
        The file name to which patient IDs should be saved. If not empty, then will save IDs of patients into another file with this name. Useful for extracting patient ID from MaCH1 output format "ID->ID". Since dataset is generally split across many files, one chromosome each, the patient IDs should be the same across these files, thus it is enough to extract the patient ID ONCE, when running this code on the smallest chromosome. For runs on all other chromosomes, leave save.ids.name="" to save time and avoid redundant work. Could name output file as "patients.fam".
}
}
\details{
        This function is needed since results of MaCH might contain weird symbols (like '2' can appear instead of A, T, C, G). This function removes all the SNPs that have not been properly imputed by MaCH, making sure that there are no missing/strange values. This is only effective when \code{\var{letter.encoding}} = True. The reason for calling this function, and not \code{\link{pre5.genos2numeric}} is because you might wish to call other software packages on the fully imputed data, which will not need the data categorized into 3 levels.  

Outputs the following files:
\preformatted{
 <file.ped>_clean<ending.ped> - in dir.out directory, the resultant file:
      the SNP columns + last columns (but no user IDs will be recorded).
 <file.dat>_clean.dat - in dir.out directory, the corresponding .dat file, will 
      be different from original <file.dat> if any bad SNPs get removed.
 <save.ids.name> - the column of patient IDs, if save.ids.name is not empty "".
}
}
\value{
        <file.ped>_clean<ending.ped> filename - the name of the output file.
}
\author{
Olia Vesselova
}

\seealso{
\code{\link{pre5.genos2numeric}}, \code{\link{pre5.genos2numeric.batch}}, 
\code{\link{pre3.call.mach}}, \code{\link{pre4.combine.case.control}}
}
\examples{
print("later")
}
