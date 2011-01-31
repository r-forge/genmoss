\name{run1.moss.regression}
\Rdversion{1.1}
\alias{run1.moss.regression}
\title{
Runs MOSS regression algorithm
}
\description{
This function performs a \acronym{MOSS} search for log-linear models as described in the paper Dobra and Massam (2009) published in Statistical Methodology.
}
\usage{
run1.moss.regression(genome.file, max.regressors = 1, chain.iterations = 10000, 
chain.replicates = 5, cutoff.max = 0.5, cutoff.min = 0.001, prob.max = 0.1, 
num.confounding.vars=0)
}
\arguments{
  \item{genome.file}{
	The input file that contains the dichotomized SNPs together with the binary outcome. The binary outcome is assumed to occupy the last column in the file.
}
  \item{max.regressors}{
	The maximum number of predictors allowed to enter a regression. Should be small, at most 5 for most applications.
}
  \item{chain.iterations}{
	The number of iterations the stochastic search algorithm will run.
}
  \item{chain.replicates}{
	The number of instances of the \acronym{MOSS} algorithm that will be run. Typically 5 instances should be run, with a minimum of 3.
}
  \item{cutoff.max}{
	The maximum Bayes factor used to determine which regressions will be retained in the list of current best models.
}
  \item{cutoff.min}{
	The minimum Bayes factor.
}
  \item{prob.max}{
	The probability of pruning of the current list of models.
}
  \item{num.confounding.vars}{
        The number of variables that must always be present in the model. These variables come after the SNP data and before the last column (which denotes the response variable).
}

}
\details{
	MOSS algorithm is run several times \code{\var{chain.replicates}}, to determine if the best regressions have been identified. After the regressions are identified, the best log-linear model associated with the regressions is found. The process outputs 6 files in the same directory where \code{\var{genome.file}} is located. These files begin with the name of the \code{\file{genome.file}}, and end with .countmodels.txt, .log, .reg, .reg.model1.txt, .spaceratio.txt, and .var. The file ending with .reg is necessary for the prediction step functions \code{\link{run2.prediction.cvv}} and \code{\link{run3.prediction.train.test}}. 
}
\value{
\item{out.name}{The name of one of the output files (ending with .reg) that will be necessary for the prediction steps.}

}
\author{
Laurent Briollais, Adrian Dobra, Olga Vesselova
}
\note{
This is a very computationally expensive step. Requires LINUX's \command{wc} functionality.
}

\seealso{
\code{\link{run2.prediction.cvv}}, \code{\link{run3.prediction.train.test}}
}
\examples{

write(rbinom(200,1,0.5), file="randbinary.txt", append=FALSE, sep=" ", ncolumns=50)
outname <- run1.moss.regression("randbinary.txt")
try(system("rm randbinary.txt*"))
}

\keyword{optimize}
\keyword{models}
\keyword{regression}
