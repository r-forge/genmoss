\name{run5.brier}
\Rdversion{1.1}
\alias{run5.brier}
\title{
Compute the Brier Score.
}
\description{
Computes the Brier score. Returns the mean and standard deviation.
}
\usage{
run5.brier(filename)
}
\arguments{
  \item{filename}{
	 The name of file that was produced either by \code{\link{run2.prediction.cvv}} or \code{\link{run3.prediction.train.test}} functions.
}
}
\value{
\item{out$mean}{The mean of the results}
\item{out$std}{The standard deviation}
}
\author{
Laurent Briollais, Adrian Dobra, Olga Vesselova
}

\seealso{
\code{\link{run1.moss.regression}}, \code{\link{run2.prediction.cvv}}, 
\code{\link{run3.prediction.train.test}}, \code{\link{run4.show.prediction}}
}
\examples{
	
write(rbinom(200,1,0.5), file="randbinary.txt", append=FALSE, sep=" ", ncolumns=50)
run1.moss.regression("randbinary.txt")
run2.prediction.cvv("randbinary.txt", "randbinary.txt.shotgun.50.1.reg")

run5.brier("randbinary.txt.49.cvv.txt")

try(system("rm randbinary.txt*"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
