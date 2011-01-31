\name{pre2.remove.genos.batch}
\Rdversion{1.1}
\alias{pre2.remove.genos.batch}
\title{
       Remove genos with many empty values for all files
}
\description{
       For all specified files, remove columns (genos) that have too many missing values. This program will automatically match CASEs and CONTROLs and their corresponding .dat files based on the specifications of prefixes, keys, and endings. 
}
\usage{
pre2.remove.genos.batch(dir.dat, dir.ped = dir.dat, dir.out, 
dir.warning = dir.out, perc.snp = 10, perc.patient = 20, empty = "0/0", 
num.nonsnp.col = 5, prefix.dat, prefix.case, prefix.control, key.dat = "", 
key.case = "CASE", key.control = "CONTROL", ending.dat = ".dat", 
ending.case = ".ped", ending.control = ".ped")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.dat}{
        The directory name where all .dat files can be found.
}
  \item{dir.ped}{
        The directory name where all .ped CASE and CONTROL files can be found. Defaults to same place as \code{\var{dir.dat}}
}
  \item{dir.out}{
        The directory name to which output files should be saved.
}
  \item{dir.warning}{
        The directory name to which warnings about patients with too many missing SNPs should go. Defaults to the same place as \code{\var{dir.out}}.
}
  \item{perc.snp}{
        The percentage (0-100 percent) of maximum empty values allowed for each geno (column). All genos that have more empty values than this threshold will be removed.
}
  \item{perc.patient}{
        The percentage (0-100 percent) of empty values allowed for each patient (row). Names of all patients who end up having more empty values than this threshold will be recorded in the warnings file.
}
  \item{empty}{
        The representation of a missing SNP value in the file ("0 0", "0/0", "1/1", "N N", etc).
}
  \item{num.nonsnp.col}{
        The number of leading columns in the .ped files that do not contain SNP values. The first columns of the file represent non-SNP values (like patient ID, gender, etc). For MaCH1 input format, the \code{\var{num.nonsnp.col=5}}, for PLINK it is 6 (due to extra disease status column).
}
  \item{prefix.dat}{
        The beginning of the file name for the .dat file (up until chrom number).
}
  \item{prefix.case}{
        The beginning of the file name for the CASE pedegree file (up until chrom number).
}
  \item{prefix.control}{
        The beginning of the file name for the CONTROL pedegree file (up until chrom number).
}
  \item{key.dat}{
        Any keyword in the name of the pedegree file that distinguishes it from other files.
}
  \item{key.case}{
        Any keyword in the name of the CASE pedegree file that distinguishes it from other non-pedegree non-CASE files.
}
  \item{key.control}{
        Any keyword in the name of the CONTROL pedegree file that distinguishes it from other non-pedegree non-CONTROL files.
}
  \item{ending.dat}{
        The ending of the .dat filenames.
}
  \item{ending.case}{
        The ending of the CASE pedegree filenames.
}
  \item{ending.control}{
        The ending of the CONTROL pedegree filenames.
}
}
\details{
         Removes SNPs that contain more than \code{\var{perc.snp}} empty geno values, from all the corresponding CASE and CONTROL .ped and .dat files in directory \code{\var{dir.dat}}. If a .ped file for some chromosome is split into several files, these files will be concatenated into one file for that chromosome, in alphabetical order. Those chromosomes that have files that satisfy the (prefix, key, ending) selection criterion but do NOT have complete set of 3 files (CASE, CONTROL, and .dat), will NOT be processed. 
}
\author{
Olia Vesselova
}


\seealso{
\code{\link{pre1.plink2mach}}, \code{\link{pre1.plink2mach.batch}}, 
\code{\link{pre2.remove.genos}}, \code{\link{pre3.call.mach}}, 
\code{\link{pre3.call.mach.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
