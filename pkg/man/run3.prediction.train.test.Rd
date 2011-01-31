\name{run3.prediction.train.test}
\Rdversion{1.1}
\alias{run3.prediction.train.test}
\title{
Prediction of test data using regression models
}
\description{
Uses the regression models identified in the \acronym{MOSS} regression step (function \code{\link{run1.moss.regression}}), to perform prediction on the \code{\var{genome.test.file}}.
}
\usage{
run3.prediction.train.test(genome.train.file, genome.test.file, models.file, 
max.regressors = 1, chain.iterations = 10000)
}
\arguments{
  \item{genome.train.file}{
	The input file that contains the dichtomized SNPs together with the binary outcome. The binary outcome is assumed to occupy the last column in the file. This is the same file as would be given as input to \acronym{MOSS} regression function \code{\link{run1.moss.regression}}.
}
  \item{genome.test.file}{
	This file is used for testing. Should be of the same format as \code{\var{genome.train.file}} with the same number of variables (columns), and arbitrary sample size (rows).
}
  \item{models.file}{
	The output file from \code{\link{run1.moss.regression}} function when run on \code{\var{genome.train.file}} with same value for \code{\var{max.regressors}} parameter. The file name would end in ".reg".
}
  \item{max.regressors}{
	Should be the same as the value used for \code{\link{run1.moss.regression}} function.
}
  \item{chain.iterations}{
	The number of samples to be drawn from the mixture of regression models.
}
}
\value{
\item{name.fitted}{The name of the output file (ending with .fitted) that will be necessary for the plotting.}
}

\author{
Laurent Briollais, Adrian Dobra, Olga Vesselova
}
\note{
Requires LINUX's \command{wc} functionality.
}

\seealso{
\code{\link{run1.moss.regression}}, \code{\link{run2.prediction.cvv}}
}
\examples{

write(rbinom(200,1,0.5), file="randbinary.txt", append=FALSE, sep=" ", ncolumns=50)
name.reg <- run1.moss.regression("randbinary.txt")

name.fitted <- run3.prediction.train.test("randbinary.txt", "randbinary.txt", name.reg)

try(system("rm randbinary.txt*"))

}
\keyword{models}
\keyword{misc}
