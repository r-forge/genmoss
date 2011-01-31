\name{get.data.dims}
\Rdversion{1.1}
\alias{get.data.dims}
\title{
Obtains matrix dimensions
}
\description{
Obtains the number of rows and columns in a matrix that is stored in a text file. The entries in the file should be either space or tab delimited. No missing values. 
}
\usage{
get.data.dims(genome.file)
}
\arguments{
  \item{genome.file}{
	Name of any file that contains a matrix of values in it, separated by either spaces or tabs. 
}
}
\value{
\item{out$nrows}{Number of rows in the matrix}
\item{out$ncols}{Number of columns in the matrix}
}
\author{
Olga Vesselova
}
\note{
Uses LINUX's \command{wc} functionality. 
}

\seealso{
\code{\link{run1.moss.regression}}, \code{\link{run2.prediction.cvv}}, 
\code{\link{run3.prediction.train.test}}
}
\examples{
write(rbinom(200,1,0.5), file="randbinary.txt", append=FALSE, sep=" ", ncolumns=50)
get.data.dims("randbinary.txt")
try(system("rm randbinary.txt*"))

}
\keyword{misc}
