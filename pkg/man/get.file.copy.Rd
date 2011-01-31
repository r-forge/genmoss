\name{get.file.copy}
\Rdversion{1.1}
\alias{get.file.copy}
\title{
        Copies files from one directory to another
}
\description{
        From given directory \code{\var{dir.in}}, copies files into \code{\var{dir.out}}. Either list of file names in \code{\var{fname}}, or all files from \code{\var{dir.in}} that start from given \code{\var{prefix}} and end with \code{\var{ending}} and contain keyword \code{\var{key}}.

}
\usage{
get.file.copy(dir.in, dir.out, fname = "", prefix = "", key = "", ending = "", 
verbal = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir.in}{
        The name of directory which contains files that need to be copied.
}
  \item{dir.out}{
        The name of directory to which files should be copied.
}
  \item{fname}{
        The list of file names (should be empty if you want it to find files itself given specifications of \code{\var{prefix}}, \code{\var{key}} and \code{\var{ending}}).
}
  \item{prefix}{
        The beginning of the file names that need to be copied.
}
  \item{key}{
        Any keyword that uniquely distinguishes the files from others.
}
  \item{ending}{
        The ending of the file names that need to be copied.
}
  \item{verbal}{
        Flag whether or not to print error messages if files with \code{\var{prefix}}, \code{\var{key}} and \code{\var{ending}} could not be found. This flag only matters if \code{\var{fname}}="".
}
}
\details{
        This function can be used in two ways:

	1. Either user provides a list of filenames that need to be copied over to \code{\var{dir.out}} directory, in which case all \code{\var{prefix}}, \code{\var{key}} and \code{\var{ending}} will be ignored.


	2. Or \code{\var{fname}}="" and some of the 3 parameters \code{\var{prefix}}, \code{\var{key}} and \code{\var{ending}} are set. In which case the program will search for files in \code{\var{dir.in}} that fulfill the specifications. 

	This function is basically file.copy, only it allows to pass in a list instead of a single file, and takes input in format that is similar to all other preprocessing functions in GenMOSS. 

}
\author{
Olia Vesselova
}

\seealso{
\code{\link{pre0.dir.create}}
}
\examples{
print("See the demo 'gendemo'.")
}
