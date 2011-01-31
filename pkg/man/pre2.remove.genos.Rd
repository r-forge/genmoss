\name{pre2.remove.genos}
\Rdversion{1.1}
\alias{pre2.remove.genos}
\title{
        Remove genos with many empty values
}
\description{
	Remove columns (genos) that have too many missing values. All genos that have more than \code{\var{perc.snp}} values missing in both \code{\var{case.ped}} AND \code{\var{control.ped}} files will be removed.
}
\usage{
pre2.remove.genos(file.dat, case.ped, control.ped, dir.dat, dir.out, 
dir.warning = dir.out, perc.snp = 10, perc.patient = 20, empty = "0/0", 
num.nonsnp.col = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.dat}{
        The name of data file as required for MaCH1.
        The file should be of the format:
         \preformatted{
         M SNP1
         M SNP2

        - Space separated
        - No header
        - Column 1: consists of "M"
        - Column 2: character SNP names
}
}
  \item{case.ped}{
        The name of pedegree data file that contains CASEs in MaCH input format.
}
  \item{control.ped}{
        The name of pedegree data file that contains CONTROLs in MaCH input format.
}
  \item{dir.dat}{
        The directory name where \code{\var{file.dat}} and \code{\var{file.ped}} can be found.
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
}
\details{
        Remove columns (genos) that have too many missing values. All genos that have more than \code{\var{perc.snp}} values missing in both \code{\var{case.ped}} AND \code{\var{control.ped}} files will be removed.


All patients that have more than \code{\var{perc.patient}} values missing will have their IDs written into "warning.<case.ped>.txt" files. Output will be two clean versions of \code{\var{case.ped}} and \code{\var{control.ped}} files in \code{\var{dir.out}} directory, and optionally the warning files in \code{\var{dir.warning}} directory.


The following files will be saved after the program is run:

\preformatted{
 - <file.dat>.removed.dat - the .dat file containing only the SNPs that were not
       removed, will be placed in dir.out directory
 - <case.ped>.removed.ped - the CASE .ped file without columns that contain too 
       many missing values based on the thresholds perc.snp; in dir.out directory
 - <control.ped>.removed.ped - the CONTROL .ped file without columns that contain
       too many missing values based on the thresholds perc.snp; 
       in dir.out directory

 - warning.<case.ped>.txt - file containing warning messages about patients that 
       have too many SNPs missing (based on perc.patients) in CASE.ped file, 
       after the removal of bad SNPs.
 - warning.<control.ped>.txt - similar to warning.<case.ped>.txt, only for 
       CONTROL file.
}

}
\author{
Olia Vesselova
}

\seealso{
\code{\link{pre1.plink2mach}}, \code{\link{pre1.plink2mach.batch}}, 
\code{\link{pre2.remove.genos.batch}}, \code{\link{pre3.call.mach}}, 
\code{\link{pre3.call.mach.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
