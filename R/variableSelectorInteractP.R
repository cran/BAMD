variableSelectorInteractP <- 
function(fname, n, p, s, nsim, keep=5, prop=0.75, 
	 codaOut="CodaChain.txt", codaIndex="CodaIndex.txt",
	 missingfile="Imputed_missing_vals", SNPsubset ) {

  commSize <- mpi.comm.size()
  allArgs <- match.call()
  allArgs[[1]] <- NULL
  
  #mpi.bcast.Robj2slave(sub1)
  #formals(variableSelector) <- allArgs
  #mpi.bcast.Robj2slave(variableSelector)
  
  A <- mpi.remote.exec(variableSelector, fname, n, p, s, nsim, keep, prop,
                       codaOut, codaIndex, missingfile, SNPsubset)
  s <- sum(SNPsubset)

  cur.table <- NULL
  for(i in 1:(commSize-1)){
    cur.table <- rbind(cur.table, A[[i]])
  }

  cur.table <- unique(cur.table)
  cur.table <- cur.table[order(cur.table[,s+1], decreasing=TRUE),]
  cur.table <- cur.table[1:keep,]

  cur.table
}
