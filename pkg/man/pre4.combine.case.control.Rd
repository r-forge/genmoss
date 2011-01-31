\name{pre4.combine.case.control}
\Rdversion{1.1}
\alias{pre4.combine.case.control}
\title{
        Combine CASE and CONTROL files
}
\description{
        Combines CASE and CONTROL files into one file, and appends disease status as the last column. The disease status is encoded as 1 for CASE and 0 for CONTROL.
}
\usage{
pre4.combine.case.control(case.file, control.file, dir.file, name.out, 
dir.out = dir.file, separ = " ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{case.file}{
        The name of the CASE file.
}
  \item{control.file}{
        The name of the CONTROL file.
}
  \item{dir.file}{
        The name of directory where CASE and CONTROL input files can be found.
}
  \item{name.out}{
        The desired name for the output file.
}
  \item{dir.out}{
        The name of directory to which output file should be written.
}
  \item{separ}{
        The separator used in the CASE and CONTROL input files.
}
}
\details{
        The function combines CASE and CONTROL together, attaching disease status as the last column: 1 for CASE and 0 for CONTROL. There will be two output files:
\preformatted{
 - <dir.out>/<name.out> - the file containing both CASE and CONTROL values, 
     with the disease status as the last column.
 - <dir.out>/<all.dat> - also will copy over ALL the files ending with ".dat" 
     that exist in \code{\var{dir.file}}.
}
}
\author{
Olia Vesselova
}


\seealso{
\code{\link{pre3.call.mach}}, \code{\link{pre3.call.mach.batch}}, 
\code{\link{pre4.combine.case.control.batch}}, 
\code{\link{pre5.genos2numeric}}, \code{\link{pre5.genos2numeric.batch}}

}
\examples{
print("See the demo 'gendemo'.")
}
