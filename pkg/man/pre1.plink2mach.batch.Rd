\name{pre1.plink2mach.batch}
\Rdversion{1.1}
\alias{pre1.plink2mach.batch}
\title{
        Convert Plink to MaCH input format for all files
}
\description{
        For all files in \code{\var{dir.in}} directory that end with \code{\var{ending.ped}} and contain keyword \code{\var{key.ped}}, and all files ending with \code{\var{ending.map}} and contain keyword \code{\var{key.map}}, runs the converter \code{\link{pre1.plink2mach}}. This will re-format all the files from Plink format to MaCH input format. 
}
\usage{
pre1.plink2mach.batch(dir.in, dir.out, ending.ped = ".ped", ending.map = ".map", 
key.ped = "", key.map = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.in}{
        The name of directory where all the files with \code{\var{ending.ped}}, \code{\var{key.ped}}, \code{\var{ending.map}} and \code{\var{key.map}} specifications are located. 
}
  \item{dir.out}{
        The name of directory to which output files .ped and .dat should be saved to.
}
  \item{ending.ped}{
        The ending of the filenames that contain pedegree data in Plink format. See format of .ped file in \code{\link{pre1.plink2mach}}.
}
  \item{ending.map}{
        The ending of the filenames that contain SNP ID information. See format of .map file in \code{\link{pre1.plink2mach}}.
}
  \item{key.ped}{
        Any keyword in the Plink pedegree file names to uniquely distinguish them from other files in the \code{\var{dir.in}} directory.
}
  \item{key.map}{
        Any keyword in the .map file names to uniquely distinguish them from other files in the \code{\var{dir.in}} directory.
}
}
\details{
        The input file formats of .ped and .map files are described in \code{\link{pre1.plink2mach}}. 
}
\author{
Olia Vesselova
}


\seealso{
\code{\link{pre1.plink2mach}}, \code{\link{pre0.dir.create}}, 
\code{\link{pre2.remove.genos}}, \code{\link{pre2.remove.genos.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
