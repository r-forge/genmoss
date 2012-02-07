\name{pre0.dir.create}
\Rdversion{1.1}
\alias{pre0.dir.create}
\title{
        Generate working subdirectory structure
}
\description{
        Function to help create the recommended subdirectory structure for the pre-processing. In dir.out a directory with name \code{\var{out.name}} will be created. Inside of this \code{\var{out.name}} directory will be a set of subdirectories, whose names will begin with \code{\var{prefix.dir}}, followed by a number, followed by short description of what the folder is designed to contain.
}
\usage{
pre0.dir.create(dir.out = ".", out.name = "newdata", prefix.dir = "d")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.out}{
        The name of directory to which new folder \code{\var{out.name}} should be saved.
}
  \item{out.name}{
        The name of the new working directory.
}
  \item{prefix.dir}{
        The start of the name of all subdirectories that will be located inside \code{\var{out.name}} folder.
}
}
\details{
        The subdirectory structure is designed to easily work with preprocessing functions of GenMOSS. Since GenMOSS preprocessing steps need to be performed in a fixed order, and there are several files per chromosome at each step, very good organization of these files is necessary to know what files have come from where and which .dat, .ped, and .fam files correspond. This function creates the directory and subdirectory structure, and it also returns the names of all the subdirectories, which can be easily used as out$d0 to out$d11. See the demo "gendemo" that shows how to effortlessly use this return variable when calling all the pre-processing steps.  
}
\value{
\item{out$d0}{The name of subdirectory into which original data should be placed.}
\item{out$d1}{The name of subdirectory into which data converted into Plink format should go. This can be done by function similar to \code{\link{ex2plink}}. }
\item{out$d2}{The name of subdirectory into which data converted into MaCH input format should go. This can be done by \code{\link{pre1.plink2mach.batch}}. }
\item{out$d3}{The name of subdirectory into which data with removed empty SNPs should go. This can be done by \code{\link{pre2.remove.genos.batch}}. }
\item{out$d4}{The name of subdirectory into which reference files needed for MaCH1 can be downloaded. }
\item{out$d5}{The name of subdirectory into which output of MaCH1 should go. This can be done by \code{\link{pre3.call.mach.batch}}. }
\item{out$d6}{The name of subdirectory into which combined CASE and CONTROL files should go. This can be done by \code{\link{pre4.combine.case.control.batch}}. }
\item{out$d7}{The name of subdirectory into which data converted to numeric 3 levels should go. This can be done by \code{\link{pre5.genos2numeric.batch}}. }
\item{out$d8}{The name of subdirectory into which data converted to binary should go. This can be done by \code{\link{pre6.discretize.batch}}. }
\item{out$d9}{The name of subdirectory into which binary data merged across all chromosomes should go. This can be done by \code{\link{pre7.merge.genos}}. }
\item{out$d10}{The name of subdirectory into which merged data split into train and test sets should go. Also all the main GenMOSS computation would go into this subdirectory. The train-test split can be done by \code{\link{pre9.split.train.test.batch}}. }
\item{out$d11}{The name of subdirectory into which desired subsets of the data should go. Also all the main GenMOSS computation for each subset would go into this subdirectory. The train-test split can be done by \code{\link{pre9.split.train.test.batch}}. }
}
\author{
Olia Vesselova
}

\seealso{
\code{\link{ex2plink}}, \code{\link{pre1.plink2mach.batch}}, 
\code{\link{pre2.remove.genos.batch}}, \code{\link{pre3.call.mach.batch}}, 
\code{\link{pre4.combine.case.control.batch}}, 
\code{\link{pre5.genos2numeric.batch}}, \code{\link{pre6.discretize.batch}}, 
\code{\link{pre7.merge.genos}}, \code{\link{pre9.split.train.test.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
