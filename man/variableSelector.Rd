\name{variableSelector}
\alias{variableSelector}
\title{ Variable Selection in Bayesian Association Model }
\description{
This function carries out variable selection on the following linear mixed model 
\deqn{ Y = X \beta + Z \gamma + \epsilon }
where the covariates 
for the random effects (in the Z-matrix) have missing values.
The Z-matrix consists of Single Nucelotide Polymorphism (SNP) data and the 
Y-vector contains the phenotypic trait of interest. The X-matrix typically 
describes the family structure of the organisms.

The best models are determined by their Bayes Factor, and uses the imputed 
values from the \code{\link{gibbsSampler}} function.
}
\usage{
variableSelector(fname, n, p, s, nsim, keep = 5, prop = 0.75,  
	codaOut="CodaChain.txt", codaIndex="CodaIndex.txt",
	missingfile = "Imputed_missing_vals", SNPsubset)
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
}
\details{
  A Metropolis-Hastings algorithm is used to conduct a stochastic search through the model space
  to find the best models.
}
\value{
  A matrix consisting of the best \code{keep} models and their Bayes Factors is returned.
}
\references{ 
  Gopal, V. "BAMD User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/assoc_model/assoc.html}
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}

Maintainer: Vik Gopal <viknesh@stat.ufl.edu>
}
\seealso{
 \code{\link{gibbsSampler}}
}
\examples{
# Load example matrices and write to csv files.
data(Y, X, Z, R, Zprob)
write.csv(cbind(Y,X,Z,R), file="generatedData.csv", quote=FALSE, row.names=FALSE)
write.csv(Zprob, file="Zprob.csv", quote=FALSE, row.names=FALSE)
	
# Run the gibbs sampler with 100 iterations, keeping the last 800
gibbsSampler(fname="generatedData.csv", fprob="Zprob.csv", n=8, p=3, s=5, nsim=1000, keep=800)

# Imputed values from gibbs sampler will be used in Variable Selector
variableSelector(fname="generatedData.csv", n=8, p=3, s=5, nsim=100, keep = 5)

#remove all generated csv files
unlink("*.csv")
unlink("*.txt")
}

\keyword{ models }
