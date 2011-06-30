\name{BAMD-package}
\alias{BAMD-package}
\alias{BAMD}
\docType{package}
\title{
Bayesian Association Model for Genomic Data with Missing Covariates
}
\description{
This package fits the following linear mixed model 
\deqn{ Y = X \beta + Z \gamma + \epsilon }
where the covariates 
for the random effects (in the Z-matrix) have missing values. 
The Z-matrix consists of Single Nucelotide Polymorphism (SNP) data and the 
Y-vector contains the phenotypic trait of interest. The X-matrix typically 
describes the family structure of the organisms.
}
\details{
\tabular{ll}{
Package: \tab BAMD\cr
Type: \tab Package\cr
Version: \tab 3.5\cr
Date: \tab 2011-06-30\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
There are two functions in this package. The first, \code{gibbsSampler}, will fit and 
estimate the posterior parameters in the model above. This will allow the 
experimenter to pick out which covariates were significantly non-zero, since 
the routine will return \eqn{(1-\alpha)100\%} confidence intervals. The imputed 
missing values at each iteration of the Gibbs sampler will be stored in a file for 
use by the second function.

The second function, \code{variableSelector}, is a variable selector that will pick out the ``best''
model, as measured by the Bayes Factor, using a stochastic search algorithm.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}

Maintainer: Vik Gopal <viknesh@stat.ufl.edu>
}
\references{
  Gopal, V. "BAMD User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/assoc_model/assoc.html}
}
\keyword{ package }
\seealso{
 \code{\link{gibbsSampler}}, \code{\link{variableSelector}}
}
