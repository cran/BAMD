# Summarize the values that were imputed for the genotypes 
# during the Gibbs sampler run.
imputedGenotypes <- function(fname, n, p, s, 
                             missingfile="Imputed_missing_vals")  {

# Extract index of missing values from input file
  DATA <- read.csv(fname, header=T)
  Y <- as.matrix(DATA[,1])
  X <- as.matrix(DATA[,2:(p+1)])
  Z <- as.matrix(DATA[,(p+2):(s+p+1)])
  rm(DATA, Y, X)
  
  missing.index <- which(Z==0, arr.ind=TRUE)
  num.missing <- dim(missing.index)[1]

  if (num.missing>0) {
    # check how many iterations of gibbs sampler were run.
    # The use of 'seek' function because of this spectacular statement in the 
    # help file for that function:
    #     "Use of ‘seek’ on Windows is discouraged.  We have found so many
    #     errors in the Windows implementation of file positioning that
    #     users are advised to use it only at their own risk, and asked not
    #     to waste the R developers' time with bug reports on Windows'
    #     deficiencies."
    dum1 <- file(missingfile, "rt")
    x <- readLines(dum1, n=-1)
    nsim.gs <- length(x)
    rm(x)
    close(dum1)

    imputedGenotypes <- matrix(0, nrow=nrow(missing.index), ncol=4)
    colnames(imputedGenotypes) <- c("recHM", "HT", "domHM", "numImputns")

    # Open file connection to imputed data
    # Compute "average" of imputed missing values
    imputed.data <- file(missingfile, "r")
    initialVals <- strsplit(readLines(imputed.data, n=1), " ")[[1]]
    imputedGenotypes[ , "numImputns"] <- 1 
    for (i in 1:nrow(missing.index)) {
     switch(initialVals[i], 
     "-1.0" = imputedGenotypes[i, "recHM"] <- imputedGenotypes[i, "recHM"] + 1,
     "0.0" = imputedGenotypes[i, "HT"] <- imputedGenotypes[i, "HT"] + 1,
     "1.0" = imputedGenotypes[i, "domHM"] <- imputedGenotypes[i, "domHM"] + 1)
    }
    for (i in 2:nsim.gs) {
      tmp <- strsplit(readLines(imputed.data, n=1), " ")[[1]]
      which.col.now <- as.integer(tmp[1])
      start.j <- min(which(missing.index[,"col"] == which.col.now))
      tmp <- tmp[-1]
      for (j in 1:length(tmp)){
       switch(tmp[j], 
         "-1.0" = imputedGenotypes[j+start.j-1, "recHM"] 
                  <- imputedGenotypes[j+start.j-1, "recHM"] + 1,
         "0.0" = imputedGenotypes[j+start.j-1, "HT"] 
                  <- imputedGenotypes[j+start.j-1, "HT"] + 1,
         "1.0" = imputedGenotypes[j+start.j-1, "domHM"] 
                  <- imputedGenotypes[j+start.j-1, "domHM"] + 1)
       imputedGenotypes[j+start.j-1, "numImputns"] <- 
         imputedGenotypes[j+start.j-1, "numImputns"] + 1
      }
    }
  
    # Close connection to imputed data
    close(imputed.data)
  
  } 
  else {
    cat("No missing values in Z-matrix!\n")
    return(NULL)
  }
  tmp.out <- imputedGenotypes[,1:3]/imputedGenotypes[,4]
  return(cbind(missing.index, tmp.out))


}
