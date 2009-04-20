vs <-
function(fname, n, p, s, nsim, keep=5, prop=0.75, betafile="beta.csv", 
	gammafile="gamma.csv", phi2file="phi2.csv", 
	sig2file="sig2.csv", missingfile="Imputed_missing_vals")  {

################################################################
# Function to update the table of models                       #
################################################################
update.table <- function(cur.table, delta, BF) {
  keep <- dim(cur.table)[1]
  BFcol <- dim(cur.table)[2]
  min.cur.BF <- min(cur.table[,BFcol])
  which.min.cur.BF <- which.min(cur.table[,BFcol])

# If updating is needed
  if(BF>min.cur.BF) {
    matching.BF <- which(BF == cur.table[,BFcol])
    if (length(matching.BF)==0) {
      cur.table[which.min.cur.BF,] <- c(delta, BF)
    }
    else{
      dup <- FALSE 
      for (i in matching.BF) {
        if (sum((cur.table[i,-BFcol] - delta)^2)==0) {
          dup <- TRUE 
        }
      }
      if (dup==FALSE) {
        cur.table[which.min.cur.BF,] <- c(delta, BF)
      }
    }
  }

  cur.table
}
################################################################



################################################################
# Function to initialise table of models                       #
################################################################
init.table <- function(keep, s) {
  cur.table <- matrix(0, nrow=keep, ncol=(s+1))
  while (length(unique(cur.table[,s+1]))<keep) {
  for (i in 1:keep ) {
    tmp <- round(runif(s))
    while((sum(tmp)==s)||(sum(tmp)==0)) {
        tmp <- round(runif(s))
    }
    cur.table[i, 1:s] <- tmp
    cur.table[i,s+1] <- computeBF(cur.table[i, 1:s])
  }
  }
  cur.table
}
################################################################


################################################################
# Function to compute estimated Bayes Factor                   #
################################################################
computeBF <- function(delta) {
	pos_delta <- which(delta>0)
	s_delta <- length(pos_delta)

	pos_L <- which(delta==0)
	s_L <- length(pos_L)

	Z_L <- Z[,pos_L]
	Z_delta <- Z[, pos_delta]
	tmp <- t(Z_L) %*% Z_L
	sqrt.det.Z_L <- sqrt(det(tmp))
	P_L <- Z_L %*% solve(tmp) %*% t(Z_L)

	rm(tmp, Z_L)

	BF <- vector("numeric", length=nsim.gs)
	sum.delta <- sum(delta)

	mean(.C("computeLoop", as.double(Y), as.double(X), as.double(Z_delta), as.double(P_L), as.double(BF), as.double(sqrt.det.Z_L), as.double(beta), as.double(gamma), as.double(phi2), as.double(sig2), as.integer(n), as.integer(p), as.integer(s), as.integer(s_delta), as.integer(pos_delta), as.integer(s_L), as.integer(pos_L), as.integer(nsim.gs), as.integer(sum.delta), PACKAGE="BAMD")[[5]])
}
################################################################

################################################################
# Function to generate a candidate delta		       #
################################################################
genDelta <- function(delta, prop) {
	tmp1 <- runif(1)
	if(tmp1<prop) {
		index <- sample(s, size=1)
		delta[index] <- 1 - delta[index]
	}
	else {
		delta <- round(runif(s))
	}
	delta
}
################################################################


# Read in X, Y, Z and R matrices
DATA <- read.csv(fname, header=T)
Y <- as.matrix(DATA[,1])
X <- as.matrix(DATA[,2:(p+1)])
Z <- as.matrix(DATA[,(p+2):(s+p+1)])
R <- as.matrix(DATA[,(s+p+2):(n+s+p+1)])
rm(DATA)


# Find the indices in Z with missing values,
# Not C-style indexing here!
missing.index <- which(Z==0, arr.ind=TRUE)
num.missing <- dim(missing.index)[1]

# read in simulated beta, gamma, phi2, sig2
beta <- t(as.matrix(read.csv(betafile, header=T, sep=",")))
gamma <- t(as.matrix(read.csv(gammafile, header=T, sep=",")))
phi2 <- as.vector(as.matrix(read.csv(phi2file, header=T, sep=",")))
sig2 <- as.vector(as.matrix(read.csv(sig2file, header=T, sep=",")))

# Extract number of simulations kept
nsim.gs <- length(sig2)

# Reparametrize to [-1 0 1]
Z[Z==1] <- -1
Z[Z==2] <-  0
Z[Z==3] <-  1

# Compute inverse sqrt of R matrix
D <- diag(eigen(R)$values)
P <- eigen(R)$vect
sqrt.R <- P %*% sqrt(D) %*% t(P)
sqrt.R.inv <- solve(sqrt.R)

if (num.missing>0) {
  # Open file connection to imputed data
  # Compute "average" of imputed missing values
  imputed.data <- file(missingfile, "r")
  ave.missing <- as.integer(strsplit(readLines(imputed.data, n=1), " ")[[1]])
  last.missing <- ave.missing
  for (i in 2:nsim.gs) {
	  tmp <- as.integer(strsplit(readLines(imputed.data, n=1), " ")[[1]]) 
	  index <- which(missing.index[,2]==tmp[1])
	  ave.missing[-index] <- ave.missing[-index] + last.missing[-index]
	  ave.missing[index] <- ave.missing[index] + tmp[-1]
	  last.missing[index] <- tmp[-1]
  }

  # Close connection to imputed data
  close(imputed.data)
  ave.missing <- ave.missing/nsim.gs

  # Find "Average" Z
  Z[missing.index] <- ave.missing
}

# Compute Ystar, Xstar and Zstar
Y <- sqrt.R.inv %*% Y
X <- sqrt.R.inv %*% X
Z <- sqrt.R.inv %*% Z

cur.table <- init.table(keep, s)

cur.delta <- rep(1, s) 
cur.delta[sample(s, size=1)] <- 0
cur.BF <- computeBF(cur.delta)

for (i in 1:nsim) {
	cand.delta <- genDelta(cur.delta, prop)
	while((sum(cand.delta)==s)||(sum(cand.delta)==0)) {
	cand.delta <- genDelta(cur.delta, prop)
	}
	cand.BF <- computeBF(cand.delta)

	MHratio <- cand.BF/cur.BF
	tmp <- runif(1)
	if(tmp<MHratio) {
		cur.delta <- cand.delta
		cur.BF <- cand.BF
                cur.table <- update.table(cur.table, cur.delta, cur.BF)
	}
	cat("Iteration ", i, " completed.\n")
}

cur.table[order(cur.table[,s+1], decreasing=TRUE),]

}

