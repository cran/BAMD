\name{variableSelectorBatchP}
\alias{variableSelectorBatchP}
\title{ Variable Selection in Parallel Batch Mode BAMD }
\description{
This function runs \code{\link{variableSelector}} in parallel 
in Batch mode.
}
\usage{
variableSelectorBatchP(fname, n, p, s, nsim, keep = 5, prop = 0.75, 
codaOut = "CodaChain.txt", codaIndex = "CodaIndex.txt", 
missingfile = "Imputed_missing_vals", SNPsubset, prefix, 
pathToLog, outfile = "out1.rdt")
}
\arguments{
  \item{fname}{ \code{fname} should be the name of a \code{.csv} file. This file should
contain the Y, X, Z and R matrices for the model, in that particular order. Hence it 
should contain \eqn{n \times (1 + p + s + n)} values. There should be a header rown in the 
input file as well. The Z matrix should use the values 1,2,3 for the SNPs and 0 for any missing SNPs. 
The program will convert the SNP codings to -1,0,1 and work with those.}
  \item{n}{ \code{n} refers to the length of the Y-vector; equivalent to the number of 
observations in the dataset. }
  \item{p}{ \code{p} is the number of columns of the X-matrix. }
  \item{s}{ \code{s} is the number of columns of the Z-matrix. Note that this is the total number of original SNPs put through the Gibbs sampler.}
  \item{nsim}{ \code{nsim} specifies the number of iterations of the Metropolis-Hastings
chain to carry out. }
  \item{keep}{ \code{keep} specifies the number of models to store. The top 
\code{keep} models will be retained. }
  \item{prop}{ As the candidate distribution for the Metropolis-Hastings chain is a mixture, one 
of whose components is a random walk, \code{prop} will determine the percentage of time that 
the random walk distribution is chosen. }
  \item{codaOut}{ This is the name of the file that was output from \code{\link{gibbsSampler}}. It contains the values obtained from the Gibbs sampler. }
  \item{codaIndex}{ This is the name of the file that describes the format of the variables in \code{codaOut}. }
  \item{missingfile}{ Contains the missing SNP values that were output from \code{\link{gibbsSampler}}. }
  \item{SNPsubset}{ A 0-1 vector of length \code{s}, indicating the SNPs that should be considered as possible variables. }
  \item{prefix}{ A prefix to name the log files from each processor. For example, if 
  prefix is specified as "rank" and there are 3 processors, then there will be 3 files 
  with names "rank00.log", "rank01.log" and "rank02.log"}
  \item{pathToLog}{ A path to where the log files should be stored. }
  \item{outfile}{ A character string - the file name to store the output table to.  }
}
\details{
  A Metropolis-Hastings algorithm is used to conduct a stochastic search through the model space
  to find the best models. \code{nsim} steps of the chain will be run on each available 
  processor. Each of them will return the best \code{keep} models they found to the 
  master. The master will strip away the duplicates and return the top 
  \code{keep} models found.

  See the scripts in demo/ directory for full examples.
}
\value{
No value is returned as it is run in Batch mode. The output object is stored in the binary output file.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\note{
Remember to copy the appropriate Rprofile that is provided in the inst/ directory to the 
directory to you are working in!
}
\seealso{
\code{\link{variableSelector}}, \code{\link{variableSelectorInteractP}}
}
\keyword{ htest }
