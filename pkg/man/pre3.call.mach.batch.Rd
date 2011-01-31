\name{pre3.call.mach.batch}
\Rdversion{1.1}
\alias{pre3.call.mach.batch}
\title{
        Call MaCH imputation with and without Hapmap
}
\description{
        This is the same program as \code{\link{pre3.call.mach}}, only it provides an easier way to set function input parameters. This is the only .batch function that does NOT run on all files. Since MaCH computation on each chromosome takes too long, it is faster to process chromosomes in parallel, rather than sequentially. This function imputes all missing values, for details, see \code{\link{pre3.call.mach}}. NOTE: In this implementation, do NOT run "with Hapmap" - so do NOT provide phases and legend files.
}
\usage{
pre3.call.mach.batch(dir.file, dir.ref = "", dir.out, prefix.dat, prefix.ped, 
prefix.phase = "", prefix.legend = prefix.phase, prefix.out = "result", 
key.dat = "", key.ped = "", key.phase = "", key.legend = "", ending.dat = ".dat", 
ending.ped = ".ped", ending.phase = ".phase", ending.legend = "legend.txt", 
chrom.num, num.iters = 2, num.subjects = 200, step2.subjects = 50, empty = "0/0", 
resample = FALSE, mach.loc = "/software/mach1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.file}{
        The name of directory where .ped and .dat files can be found. The format of these files is described in \code{\link{pre3.call.mach}}
}
  \item{dir.ref}{
        The name of directory where .phase and .legend files have been downloaded to.
}
  \item{dir.out}{
        The name of directory to which output files should go.
}
  \item{prefix.dat}{
        The beginning of the file name for the .dat file (up until chrom number).
}
  \item{prefix.ped}{
        The beginning of the file name for the .ped pedegree file (up until chrom number).
}
  \item{prefix.phase}{
        The beginning of the file name for the phase file (up until chrom number). This file can be obtained from websites like: \url{http://hapmap.ncbi.nlm.nih.gov/downloads/phasing/2007-08_rel22/phased/} or similar/updated versions. No zip. Must be a normal and readable by R file.
}
  \item{prefix.legend}{
        The beginning of the file name for the legend file (up until chrom number). This file can be obtained from same website as phase file. No zip.
}
  \item{prefix.out}{
        The prefix for naming output files that MaCH1 should use. If \code{\var{num.subjects}} > 0 then the \code{\var{num.subjects}} will be appended to the prefix name.
}
  \item{key.dat}{
        Any keyword in the name of the .dat file that distinguishes it from other files.
}
  \item{key.ped}{
        Any keyword in the name of the pedegree file that distinguishes it from other files.
}
  \item{key.phase}{
        Any keyword in the name of the phase file that distinguishes it from other files.
}
  \item{key.legend}{
        Any keyword in the name of the legend file that distinguishes it from other files.
}
  \item{ending.dat}{
        The ending of the .dat filename.
}
  \item{ending.ped}{
        The ending of the pedegree filename.
}
  \item{ending.phase}{
        The ending of the phase filename.
}
  \item{ending.legend}{
        The ending of the legend filename.
}
  \item{chrom.num}{
        The chromosome number for which processing should be done.
}
  \item{num.iters}{
        The how many iterations MaCH1 should make in its first step to estimate its model parameters.
}
  \item{num.subjects}{
        How many individuals from the sample should be used for model building by the first step of MaCH1. The random subset of inidividuals will be extracted by this program. Recommended number of subjects is 200-500. Value <= 0 corresponds to using ALL the subjects in the dataset.
}
  \item{step2.subjects}{
        How many individuals should be processed at a time during the second step of MaCH computation. Value <= 0 will use ALL the subjects in the dataset. This variable is important to reduce exponential computation time required by MaCH when number of individuals is too large.

}
  \item{empty}{
        The way a missing/empty entry of SNP is represented in pedegree file.
}
  \item{resample}{
        Whether or not to overwrite the existing file containing the \code{\var{num.subjects}} entries produced by previous runs of this algorithm with same .dat, .ped, and \code{\var{num.subjects}} parameters. By default, if the subjects have been sampled before, they are re-used.

}
  \item{mach.loc}{
        The location directory where "mach" executable can be found.
}
}
\details{
        This function imputes all missing values in the data. See \code{\link{pre3.call.mach}} for details. This is the same program as \code{\link{pre3.call.mach}}, only it provides an easier way to set function input parameters. Recall that \code{\link{pre3.call.mach}} function requires users to specify names of .ped, .dat, .phase, and .legend for each chromosome - these files normally would have exactly same names across all chromosomes, and would only differ by the chromosome number. Thus after running \code{\link{pre3.call.mach}}, for chromosome 1, and in order to run next chromosome (say, chrom "2"), user would need to change this chromosome number in 4 places: from "1" to "2" in .ped, .dat, .phase, and .legend. This function allows user to just change one variable \code{\var{chrom.num}}, from "1" to "2", and all the other files will be obtained automatically.


 This is the only .batch function that does NOT run on all files. Since MaCH computation on each chromosome takes too long, it is faster to process chromosomes in parallel, rather than sequentially. Thus if your dataset is large, then it is recommended to run this function on different computers/nodes for different chromosomes.
}
\references{
       MaCH website: \url{http://www.sph.umich.edu/csg/abecasis/MACH/download/} 
}
\author{
Olia Vesselova
}
\note{
       In this current version, avoid using Hapmap. So do NOT provide reference and legend files. 
}


\seealso{
\code{\link{pre2.remove.genos}}, \code{\link{pre2.remove.genos.batch}}, 
\code{\link{pre3.call.mach}}, \code{\link{pre4.combine.case.control}}, 
\code{\link{pre4.combine.case.control.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
