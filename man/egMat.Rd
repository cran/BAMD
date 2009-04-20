\name{egMat}
\alias{Y}
\alias{X}
\alias{Z}
\alias{R}
\alias{Zprob}
\docType{data}
\title{ Example matrices for input to gs() and vs()}
\description{
  These are example datasets that demonstrate the format of input 
  to the routines in BAMD.
}
\usage{data(Y)}
\format{
  For this dataset, \code{n}=8, \code{p}=3 and \code{s}=5.
}
\references{
  Gopal, V. "BAMD User Manual"
  \url{http://www.stat.ufl.edu/~viknesh/assoc_model/assoc.html}
}
\examples{
data(Y,X,Z,R,Zprob)
write.csv(cbind(Y,X,Z,R), file="generatedData.csv", quote=FALSE, row.names=FALSE)
write.csv(Zprob, file="Zprob.csv", quote=FALSE, row.names=FALSE)
}
\keyword{datasets}
