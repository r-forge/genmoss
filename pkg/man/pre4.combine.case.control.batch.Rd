\name{pre4.combine.case.control.batch}
\Rdversion{1.1}
\alias{pre4.combine.case.control.batch}
\title{
        Combine CASE and CONTROL files for all files
}
\description{
        For each pair of CASE and CONTROL files, combine them into one file. Last column of each output file will contain the disease status. The disease status is encoded as 1 for CASE and 0 for CONTROL.
}
\usage{
pre4.combine.case.control.batch(dir.file, dir.out = dir.file, prefix.case, 
prefix.control, prefix.out, key.case = "", key.control = "", 
ending.case = ".mlgeno", ending.control = ".mlgeno", separ = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.file}{
        The name of directory where CASE and CONTROL files can be found.
}
  \item{dir.out}{
        The name of directory to which output file should be written.
}
  \item{prefix.case}{
        The beginning of the file name for the CASE file (up until chrom number).
}
  \item{prefix.control}{
        The beginning of the file name for the CONTROL file (up until chrom number).
}
  \item{prefix.out}{
        The beginning of the file name for the output file (up until chrom number).
}
  \item{key.case}{
        Any keyword in the name of the CASE file that distinguishes it from other files.
}
  \item{key.control}{
        Any keyword in the name of the CONTROL file that distinguishes it from other files.
}
  \item{ending.case}{
        The ending of the CASE filename.
}
  \item{ending.control}{
        The ending of the CONTROL filename.
}
  \item{separ}{
        The separator used in the CASE and CONTROL input files.
}
}
\details{
        The function combines CASE and CONTROL together, attaching disease status as the last column: 1 for CASE and 0 for CONTROL. There will be two output files for each pair of CASE and CONTROL:
\preformatted{
 - <dir.out>/<prefix.out><chrom.num><ending.case> - the file containing both 
     CASE and CONTROL values, with the disease status as the last column.
 - <dir.out>/<all.dat> - also will copy over ALL the files ending with ".dat" 
     that exist in \code{\var{dir.file}}.
}

}
\author{
Olia Vesselova
}


\seealso{
\code{\link{pre3.call.mach}}, \code{\link{pre3.call.mach.batch}}, 
\code{\link{pre4.combine.case.control}}, \code{\link{pre5.genos2numeric}}, 
\code{\link{pre5.genos2numeric.batch}}
}
\examples{
print("See the demo 'gendemo'.")
}
