`gibbsSampler` <-
function(fname, fprob, n, p, s, a=2.01, b=0.99099, c=2.01, d=0.99099, 
		beta.in=1, gamma.in=1, sig2.in=1, phi2.in=1, alpha=0.05,
		nsim=1000, keep=100, codaOut="CodaChain.txt", 
		codaIndex="CodaIndex.txt") {
DATA <- read.csv(fname, header=T)
Y <- as.matrix(DATA[,1])
X <- as.matrix(DATA[,2:(p+1)])
Z <- as.matrix(DATA[,(p+2):(s+p+1)])
R <- as.matrix(DATA[,(s+p+2):(n+s+p+1)])
rm(DATA)

if(missing(fprob)) {
# change this to 0 or your birthdate!
#  PROB <- matrix(1/3, nrow=n, ncol=3*s)
  SNPprior = 0
  PROB.sub <- 0
}
else {
  PROB <- read.csv(fprob, header=T)
  PROB <- as.matrix(PROB)
  SNPprior = 1
}

# Find the indices in Z with missing values,
# 1 is subtracted from missing.index because we want C-style indexing.
missing.index <- which(Z==0, arr.ind=TRUE)
num.missing <- dim(missing.index)[1]

  if(num.missing>0) {
    if(SNPprior==1) {
      PROB.sub <- NULL
      for (i in 1:num.missing) 
        PROB.sub <- rbind(PROB.sub, PROB[missing.index[i,1], (3*(missing.index[i,2]-1)+1) : (3*missing.index[i,2])])
      rm(PROB)
    }
    else PROB.sub <- 0
  }
  else { 
    PROB.sub <- 0
    if (SNPprior==1) rm(PROB)
  }
  missing.index <- missing.index - 1


# Re-code the 1's to -1's, 2's to 0's and 3's to -1's
# Hence automatically, the missing values are coded with SNP value "0" 
# to begin with.
# Coerce Z to a matrix because it is a list to begin with.
  Z[Z==1] <- -1
  Z[Z==2] <-  0
  Z[Z==3] <-  1

# Set up the matrices to hold the simulated data
# NOTE that the beta_s and gamma_s matrices are p by nsim and 
# s by nsim, because we need to pass a contiguous part of memory to the 
# component C functions
beta_s <- matrix(0, nrow=p, ncol=nsim)
gamma_s <- matrix(0, nrow=s, ncol=nsim)
sig2_s <- vector(mode="numeric", length=nsim)
phi2_s <- vector(mode="numeric", length=nsim)

# Initialize the above vectors/matrices and anything else needed
if(!missing(beta.in)) {
  if(length(beta.in)==p)
  beta_s[,1] <- beta.in
  else 
  stop("Length of initial beta vector is not equal to 'p'. Please check!")
}
else
  beta_s[,1] <- 1

if(!missing(gamma.in)) {
  if(length(gamma.in)==s)
  gamma_s[,1] <- gamma.in
  else 
  stop("Length of initial gamma vector is not equal to 's'. Please check!")
}
else
  gamma_s[,1] <- 1

if(!missing(sig2.in)) {
  if(length(sig2.in)==1)
  sig2_s[1] <- sig2.in
  else 
  stop("Initial sig2 value is not a scalar. Please check!")
}
else
  sig2_s[1] <- 1

if(!missing(phi2.in)) {
  if(length(phi2.in)==1)
  phi2_s[1] <- phi2.in
  else 
  stop("Initial phi2 value is not a scalar. Please check!")
}
else
  phi2_s[1] <- 1

#dyn.load("gibbs.so")

param <- 1
cat("Loading compiled code..\n")
tmp <- .C("gibbs_sampler",  as.double(Y), as.double(X), as.double(Z), as.double(R), as.double(PROB.sub),
as.integer(missing.index), as.integer(num.missing), as.integer(nsim), as.integer(n), as.integer(p), 
as.integer(s), as.double(a), as.double(b), as.double(c), as.double(d), as.double(beta_s), 
as.double(gamma_s), as.double(sig2_s), as.double(phi2_s), as.integer(param), as.integer(keep),
as.integer(SNPprior), PACKAGE="BAMD")
#dyn.unload("gibbs.so")
cat("Computations by compiled code finished!\n")

# Only keep the final "keep" values!
store <- nsim-keep+1
phi2.keep <- tmp[[19]][store:nsim]
sig2.keep <- tmp[[18]][store:nsim]
store <- (nsim-keep)*s +1
gamma.keep <- t(matrix(tmp[[17]][store:(nsim*s)],nrow=s,ncol=keep))
store <- (nsim-keep)*p +1
beta.keep <- t(matrix(tmp[[16]][store:(nsim*p)],nrow=p,ncol=keep))

# write CODA index file
zz <- file(codaIndex, "w")
for ( i in 1:p ) {
  tmp <- paste("beta","[",i,"]", sep="")
  tmp <- paste(tmp, (i-1)*keep+1, i*keep, "\n", sep=" ")
  cat(tmp, file=zz, append=TRUE)
}
for ( i in 1:s ) {
  tmp <- paste("gamma","[",i,"]", sep="")
  tmp <- paste(tmp, (p+i-1)*keep+1, (p+i)*keep, "\n", sep=" ")
  cat(tmp, file=zz, append=TRUE)
}
  tmp <- paste("phi2", (p+s)*keep+1, (p+s+1)*keep, "\n", sep=" ")
  cat(tmp, file=zz, append=TRUE)
  tmp <- paste("sig2", (p+s+1)*keep+1, (p+s+2)*keep, "\n", sep=" ")
  cat(tmp, file=zz, append=TRUE)
close(zz)

# write CODA chain file
tmp <- c(as.vector(beta.keep), as.vector(gamma.keep), phi2.keep, sig2.keep)
write.table(cbind(rep(1:keep, times=(p+s+2)), round(tmp, digits=5)), 
	file=codaOut, quote=FALSE, col.names=FALSE, row.names=FALSE)

}

