## Interactive code: Make sure slaves are already spawned!
##                   Make sure the right profile file is used!
  library(BAMD, lib.loc="/scratch/ufhpc/viknesh/Rlibs")
  mpi.remote.exec(library(BAMD, lib.loc="/scratch/ufhpc/viknesh/Rlibs"))
 
Y <- matrix(c(-66.0665,1,0,0,2,2,2,3,0,1,-0.88,0,0,0,0,0,0,
-298.342,1,0,0,3,2,2,1,0,-0.88,1,0,0,0,0,0,0,
155.0717,0,1,0,3,3,2,3,0,0,0,1,0,0,0,0,0,
442.3504,0,1,0,1,2,3,1,3,0,0,0,1,0,0,0,0,
7.7177,0,1,0,3,2,1,2,3,0,0,0,0,1.25,0,0,0,
810.5861,0,1,0,1,3,3,1,1,0,0,0,0,0,1,0,0,
49.1489,0,0,1,2,0,1,3,1,0,0,0,0,0,0,1,-0.13,
451.2934,0,0,1,0,1,3,2,3,0,0,0,0,0,0,-0.13,1), nrow=8, byrow=TRUE)
write.csv(Y, file="gen.data.20.csv", quote=FALSE, row.names=FALSE) 

 
  sub1 <- c(1,0,1,1,0)

  gibbsSampler("gen.data.20.csv", n=8, p=3, s=5, nsim=10000, keep=7000)

  A <- variableSelectorInteractP(fname="gen.data.20.csv", n=8, p=3, s=5, 
        nsim=6, keep=5, prop=0.75, codaOut="CodaChain.txt", 
        codaIndex="CodaIndex.txt", missingfile="Imputed_missing_vals",
        SNPsubset=sub1)

 unlink("gen.data.20.csv")
 mpi.close.Rslaves(dellog=FALSE)
 mpi.quit()

######## PARALLEL BATCH MODE CODE ##############
## testing parallel code
## make sure the Coda output and imputed missing values are already there!

#  library(BAMD, lib.loc="/scratch/ufhpc/viknesh/Rlibs")
#  variableSelectorBatchP("gen.data.20.csv", n=8, p=3, s=5, nsim=7, keep=3, 
#                       SNPsubset=c(0,1,1,1,0), prefix="rankTest", 
#                       pathToLog="log/", outfile="finaltest.rdt")
#  mpi.quit()
