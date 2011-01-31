\name{pre3.call.mach}
\Rdversion{1.1}
\alias{pre3.call.mach}
\title{
        Call MaCH imputation with and without Hapmap 
}
\description{
	Calls MACH1 program on \code{\var{file.ped}} and \code{\var{file.dat}}. MaCH1 can be run in 2 different ways: 1. with Hapmap, and 2. without Hapmap. NOTE: In this implementation, do NOT run "with Hapmap".

This program first runs MaCH1 on \code{\var{file.ped}} with Hapmap to fill in missing values for those SNPs that exist in the reference file; and then MaCH1 is run on the result without Hapmap to fill in all the remaining missing values. If no reference files \code{\var{ref.phase}} and \code{\var{ref.legend}} are provided, then the program runs MaCH1 without Hapmap only. To clean up any weird MaCH output, use \code{\link{genos.clean}} or \code{\link{pre5.genos2numeric}}.

}
\usage{
pre3.call.mach(file.dat, file.ped, dir.file, ref.phase = "", ref.legend = "", 
dir.ref = "", dir.out, out.prefix = "result", chrom.num = "", num.iters = 2, 
num.subjects = 200, step2.subjects = 50, empty = "0/0", resample = FALSE, 
mach.loc = "/software/mach1")
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
  \item{file.ped}{
        The name of pedegree data file in MaCH1 input format.
\preformatted{
 p1    p1    0       0       1       C/C     N/N     T/C ...
 p2    p2    0       0       1       T/T     A/C     G/G ...
 ...

   - Tab separated
   - Alleles are separated by slash '/' (IMPORTANT!)
   - No header
   - 5 non-SNP leading columns
   - Col 1: sample/patient ID: some unique ID
   - Col 2: family ID: can be same as patient ID
   - Col 3 and Col 4: parents: mother/father: can all be 0
   - Col 5: gender, 1-male, 2-female
   - Col 6+: geno information, slash separator between alleles.
}


}
  \item{dir.file}{
        The name of directory where \code{\var{file.ped}} can be found.
}
  \item{ref.phase}{
        The name of the reference file, must have no missing values, can be obtained from websites like: \url{http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2007-08_rel22/phased/} or similar/updated versions. No zip. Must be a normal and readable by R file.
}
  \item{ref.legend}{
        The name of legend file for \code{\var{file.phase}}, obtained from same website. No zip.
}
  \item{dir.ref}{
        The name of directory where \code{\var{ref.phase}} and \code{\var{ref.legend}} can be found.
}
  \item{dir.out}{
        The name of directory where MaCH1 output should go.
}
  \item{out.prefix}{
        The prefix for naming output files that MaCH1 should use. If \code{\var{num.subjects}} > 0 then the \code{\var{num.subjects}} will be appended to the prefix name.
}
  \item{chrom.num}{
        The optional string denoting the chromosome number, for better naming of intermediate files.
}
  \item{num.iters}{
        The number of iterations MaCH1 should make in its first step to estimate its model parameters. The same number will be used for parameter estimation when using Hapmap and when NO Hapmap is used. 
}
  \item{num.subjects}{
        How many individuals from the sample should be used for model building by the first step of MaCH1. The random subset of inidividuals will be extracted by this program. Recommended number of subjects is 200-500. Value <= 0 corresponds to using ALL the subjects in the dataset.
}
  \item{step2.subjects}{
        How many individuals should be processed at a time during the second step of MaCH computation. Value <= 0 will use ALL the subjects in the dataset. This variable is important to reduce exponential computation time required by MaCH when number of individuals is too large. However if this number is too low, the second step of MaCH might not get enough samples, thus making weird prediction of '2' instead of an Allele value. To reduce the number of '2's, try to set step2.subjects to a larger value. To remove all SNPs that have a 2 predicted for any of its entries, use \code{\link{genos.clean}} or \code{\link{pre5.genos2numeric}}.
}
  \item{empty}{
        The way a missing/empty entry of SNP is represented in \code{\var{file.ped}}.
}
  \item{resample}{
        Whether or not to overwrite the existing file containing the \code{\var{num.subjects}} entries produced by previous runs of this algorithm with same \code{\var{file.dat}}, \code{\var{file.ped}} and \code{\var{num.subjects}} parameters. By default, if the subjects have been sampled before, they are re-used.
}
  \item{mach.loc}{
        The location directory where "mach" executable can be found.
}
}
\details{
        This program first runs MaCH1 on \code{\var{file.ped}} with Hapmap to fill in missing values for those SNPs that exist in the reference file; and then MaCH1 is run on the result without Hapmap to fill in all the remaining missing values. If no reference files \code{\var{ref.phase}} and \code{\var{ref.legend}} are provided, then the program runs MaCH1 without Hapmap only.

It is recommended to avoid using Hapmap functionality in this implementation. 

        The MaCH1 algorithm requires 2 steps to be performed. The first step of MaCH1 will be run on num.subjects randomly chosen from the set. The file with randomly chosen individuals will be saved as \code{\var{file.ped.<num.subjects>.ped}} in \code{\var{dir.file}} directory. If the file already exists for this num.subjects, the old file will be used if \code{\var{resample}}=F. If \code{\var{resample}}=T then old files will be ignored, and new sampling will take place. The step1 of MaCH will only be run if \code{\var{resample}}=T, or if the files that MaCH1 produces do not exist yet. Thus if step1 runs well, but step2 crashes, re-calling this function will not waste time on re-running step1 over again.

The second step without Hapmap takes exponentially long wrt number of subjects processed. Thus the second step will be run on bunches of subjects, \code{\var{step2.subjects}} at a time.

A subdirectory structure for debugging will be formed in \code{\var{dir.out}}, the directory will be named 'working'. 

Two output files will be produced in \code{\var{dir.out}}: the .ped file that will not have any missing values, will be named <\code{\var{out.prefix}}><\code{\var{chrom.num}}>.mlgeno, and a .dat file (same as before). 

}
\references{
       MaCH website: \url{http://www.sph.umich.edu/csg/abecasis/MACH/download/} 
}
\author{
Olia Vesselova
}
\note{
        Since instead of filling in missing values, MaCH1 is re-predicting ALL the values in the dataset, the Hapmap functionality is not desirable. Thus avoid using Hapmap reference files. 


	Also, MaCH prediction is not always valid, as it may contain Allele of value '2' (when only A, C, T, G are used). Programs \code{\link{pre5.genos2numeric}} and \code{\link{genos.clean}} help to remove those. 
}


\seealso{
\code{\link{pre2.remove.genos}}, \code{\link{pre2.remove.genos.batch}}, 
\code{\link{pre3.call.mach.batch}}, \code{\link{pre4.combine.case.control}},
\code{\link{pre4.combine.case.control.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
