\name{run2.prediction.cvv}
\Rdversion{1.1}
\alias{run2.prediction.cvv}
\title{
Perform cross-validation on regression models
}
\description{
Uses the regression models identified in the \acronym{MOSS} regression step (function \code{\link{run1.moss.regression}}), to perform prediction by cross-validation.
}
\usage{
run2.prediction.cvv(genome.file, models.file, max.regressors = 1, cvv.fold = 2, 
chain.iterations = 10000)
}
\arguments{
  \item{genome.file}{
	The input file that contains the dichtomized SNPs together with the binary outcome. The binary outcome is assumed to occupy the last column in the file. This is the same file as would be given as input to MOSS regression function \code{\link{run1.moss.regression}}.
}
  \item{models.file}{
	The output file from \code{\link{run1.moss.regression}} function when run on \code{\var{genome.file}} with same value for \code{\var{max.regressors}} parameter. This filename ends with .reg. 
}
  \item{max.regressors}{
	Should be the same as the value used for \code{\link{run1.moss.regression}} function.
}
  \item{cvv.fold}{
	The type of cross-validation to perform. For example cvv.fold=2 splits the samples in half, fits the models with one half and predicts the second half. cvv.fold=k performs k-fold validation.
}
  \item{chain.iterations}{
	The number of samples to be drawn from the mixture of regression models.
}
}
\details{
Performs cross-validation on the \code{\var{genome.file}}, writes a file with name ending with ".cvv.txt". 
}
\value{
\item{name.cvv}{The name of the output file (ending with .cvv.txt) that will be necessary for the plotting.}
}
\author{
Laurent Briollais, Adrian Dobra, Olga Vesselova
}
\note{
Requires LINUX's \command{wc} functionality.
}

\seealso{
\code{\link{run1.moss.regression}}, \code{\link{run3.prediction.train.test}}
}
\examples{

write(rbinom(200,1,0.5), file="randbinary.txt", append=FALSE, sep=" ", ncolumns=50)
name.reg <- run1.moss.regression("randbinary.txt")

name.cvv <- run2.prediction.cvv("randbinary.txt", name.reg)

try(system("rm randbinary.txt*"))

}

\keyword{models}
\keyword{misc}
