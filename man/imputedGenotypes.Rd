\name{imputedGenotypes}
\alias{imputedGenotypes}
\title{ Estimate posterior distribution of missing SNPs }
\description{
After the Gibbs sampler has been run, this function retrieves the values 
that were imputed for the missing SNP values and prints out the frequency 
with which the particular values were sampled. Assuming the chain has 
reached stationarity (that sufficient burn-in has been discarded), this 
corresponds to the posterior distribution of the SNPs.
}
\usage{
imputedGenotypes(fname, n, p, s, missingfile = "Imputed_missing_vals")
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
  \item{s}{ \code{s} is the number of columns of the Z-matrix. }
  \item{missingfile}{ Contains the missing SNP values that were output from \code{\link{gibbsSampler}}. }
}
\details{
  This posterior distribution is different from the "average" Z-matrix that 
  is used in the variable selector method. In that situation, we ignore the 
  fact that at every iteration of the Gibbs sampler, only one column is 
  updated and the majority of the Z-matrix remains the same. 

  In this situation, we consider a sample of the genotype to be new only 
  every time it is updated, not every Gibbs sampler iteration. For example,
  if there were 2 columns with missing values and 7000 Gibbs sampler 
  iterations were kept, only ~3500 samples would be taken for each column. 
}
\value{
  A matrix specifying the posterior distribution for each missing SNP.
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
imputedGenotypes("generatedData.csv", n=8, p=3, s=5)

#remove all generated csv files
unlink("*.txt")
unlink("*.csv")
}
\keyword{ models }
