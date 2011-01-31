\name{pre1.plink2mach}
\Rdversion{1.1}
\alias{pre1.plink2mach}
\title{
        Convert Plink to MaCH input format
}
\description{
        Provided with Plink-format files \code{\var{file.ped}} and \code{\var{file.map}} in \code{\var{dir.in}}, this function re-formats it into MACH pedigree (\code{\var{file.ped}}) and data (\code{\var{file.dat}}) file formats, and saves the reformatted files in \code{\var{dir.out}}.
}
\usage{
pre1.plink2mach(file.ped = "", file.map = "", dir.in, dir.out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.ped}{
        The name of the pedegree file. This file should be in Plink format:

\preformatted{
 p1    p1    0       0       1       2       C/C     N/N     T/C ...
 p2    p2    0       0       1       2       T/T     A/C     G/G ...
 ...

   - Tab separated
   - No header
   - 6 non-SNP leading columns
   - Col 1 and Col 2: patient ID: some unique ID
   - Col 3 and Col 4: parents: mother/father: can be set to 0
   - Col 5: gender, 1 - male, and 2 - female
   - Col 6: disease status: 1 CONTROL and 2 CASE
   - Col 7+: geno information, slash separator between alleles.
}

}
  \item{file.map}{
        The name of the .map input file. This file should contain names of SNPs in the following format:
\preformatted{
 19 rs32453434 0 5465475
 19 rs6547434 0 23534543
 ...

   - Space separated
   - No header
   - 4 columns:
   - Col 1: Chromosome number
   - Col 2: SNP ID or any other marker for SNP
   - Col 3: genetic distance (can be set to 0)
   - Col 4: physical locations (can be set to 0)
   - Number of rows is the number of SNPs used in the given 
       chromosome. (= number of SNP columns of .ped)
}


}
  \item{dir.in}{
        The directory where \code{\var{file.ped}} and \code{\var{file.map}} can be found.
}
  \item{dir.out}{
        The directory to which output .ped and .dat files should go.
}
}
\details{
This function converts from Plink to MaCH input format. There is no need to specify both \code{\var{file.ped}} and \code{\var{file.map}}; so one of them can be an empty string (""), in which case, this file will not be processed. So that you can use this function to do ONLY PED files but not map, and vice versa. 

}
\author{
Olia Vesselova
}
\note{
        Note: the function does NOT change unknown Allele values from "0" to "N", as MACH program can use either. Does NOT recode gender to "M" and "F", since MaCH1 doesn't care, but further file processing interprets "F" as "FALSE".

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pre1.plink2mach.batch}}, \code{\link{pre0.dir.create}}, 
\code{\link{pre2.remove.genos}}, \code{\link{pre2.remove.genos.batch}} 
}
\examples{
print("See the demo 'gendemo'.")
}
