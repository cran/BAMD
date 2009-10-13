\name{gibbsSampler}
\alias{gibbsSampler}
\title{ Estimate posterior parameters in Bayesian Association Model }
\description{
This function fits the following linear mixed model 
\deqn{ Y = X \beta + Z \gamma + \epsilon }
where the covariates 
for the random effects (in the Z-matrix) have missing values.
The Z-matrix consists of Single Nucelotide Polymorphism (SNP) data and the 
Y-vector contains the phenotypic trait of interest. The X-matrix typically 
describes the family structure of the organisms.

The model is fit by embedding it in a Bayesian framework and estimating the 
posterior parameters using a Gibbs sampler.
}
\usage{
gibbsSampler(fname, fprob, n, p, s, a = 2.01, b = 0.99099, c = 2.01, d = 0.99099, 
	beta.in = 1, gamma.in = 1, sig2.in = 1, phi2.in = 1, alpha = 0.05, 
	nsim = 1000, keep = 100, codaOut="CodaChain.txt", 
	codaIndex="CodaIndex.txt")
}
\arguments{
  \item{fname}{ \code{fname} should be the name of a \code{.csv} file. This file should
contain the Y, X, Z and R matrices for the model, in that particular order. Hence it 
should contain \eqn{n \times (1 + p + s + n)} values. There should be a header rown in the 
input file as well. The Z matrix should use the values 1,2,3 for the SNPs and 0 for any missing SNPs. 
The program will convert the SNP codings to -1,0,1 and work with those.}
  \item{fprob}{ \code{fprob} should also be a \code{.csv} file. It should contain one 
probability vector for each entry in the Z-table. Hence it should be a matrix of 
dimension \eqn{n \times 3s}. The program will read in the entire table, but only store 
the distributions corresponding to the missing values. If uniform priors are to be used, 
there is no need to specify anything. }
  \item{n}{ \code{n} refers to the length of the Y-vector; equivalent to the number of 
observations in the dataset. }
  \item{p}{ \code{p} is the number of columns of the X-matrix. }
  \item{s}{ \code{s} is the number of columns of the Z-matrix. }
  \item{a,b,c,d}{ \code{ab,c,d} are hyperparameters in the Bayesian set-up. }
  \item{beta.in}{ \code{beta.in} is the initial value for the Gibbs sampler. It should be 
a vector of length \code{p}. }
  \item{gamma.in}{ \code{gamma.in} is the initial value for the Gibbs sampler. It should 
be a vector of length \code{s}. }
  \item{sig2.in}{ \code{sig2.in} is the initial value for the Gibbs sampler. It should be 
a vector of length 1. }
  \item{phi2.in}{ \code{sig2.in} is the initial value for the Gibbs sampler. It should be 
a vector of length 1. }
  \item{alpha}{ \code{alpha} refers to the \eqn{(1- \alpha)100\%} confidence intervals 
that the program should output. }
  \item{nsim}{ \code{nsim} specifies the number of iterations of the Gibbs sampler 
to carry out. }
  \item{keep}{ \code{keep} specifies which values from the Gibbs sampler chain to keep 
and use when computing the mean and confidence intervals. This allows user to allow 
for burn-in. }
  \item{codaOut}{ File to store Gibbs sample parameters in coda format}
  \item{codaIndex}{ File to store description of sampled and stored values. }
}
\details{
  For further details on the prior distributions used, please refer to the User Guide 
  in the reference(s) given below.
}
\value{
  There will be no \code{R} object returned. Instead, as the routine is running, it 
  will print debugging statements to show the user which iteration of the Gibbs sampler it is 
  currently at. This would allow the user to detect if something is going wrong with the routine. 

  The values sampled from the full conditionals will be stored in 
  the files \code{CodaChain.txt} and \code{Imputed_missing_vals} 
  by default. 
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
 \code{\link{variableSelector}}
}
\examples{
# Load example matrices and write to csv files.
data(Y, X, Z, R, Zprob)
write.csv(cbind(Y,X,Z,R), file="generatedData.csv", quote=FALSE, row.names=FALSE)
write.csv(Zprob, file="Zprob.csv", quote=FALSE, row.names=FALSE)
	
# Run the gibbs sampler with 100 iterations, keeping the last 800
gibbsSampler(fname="generatedData.csv", fprob="Zprob.csv", n=8, p=3, s=5, nsim=1000, keep=800)

#remove all generated csv files
unlink("*.txt")
unlink("*.csv")
}
\keyword{ models }
