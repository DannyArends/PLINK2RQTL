# Convert .ped & .map files to R/qtl csvr format
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

#Set this to where YOU stored the example data
setwd("D:/Github/PLINK2RQTL/test")

convertToCSVR <- function(ped = "test.ped", map = "test.map", out = "cross.csvr", verbose = FALSE){
  mapdata <- read.table(map, colClasses=c("character"))                                                             # Read the MAP data
  colnames(mapdata) <- c("Chr", "ID", "cM", "BP")

  SNPcolnames <- paste0(unlist(lapply(mapdata[,"ID"],rep,2)), c(".A",".B"))                                         # Create column names for the SNPs

  peddata <- read.table(ped, colClasses=c("character"), na.strings=c("-9"))                                         # Read the PED data
  colnames(peddata) <- c("FID", "IID", "PID","MID", "Sex", "Pheno", SNPcolnames)

  genotypes <- matrix(NA, length(mapdata[,"ID"]), nrow(peddata))                                                    # Empty genotype matrix
  rownames(genotypes) <- mapdata[,"ID"]
  for(snp in mapdata[,"ID"]){
    cols <- grep(snp, colnames(peddata))                                                                            # Get the columns associated with this SNP
    snpalleles <- sort(unique(unlist(peddata[,cols])))                                                              # The SNP alleles
    if(verbose) cat("For", snp,"found", snpalleles,"\n")                                                            # Some info / debug
    genotype <- apply(peddata[,cols],1,function(x){
      # NOTE: I do NOT handle missing genotype data (YET)
      if(x[1] != x[2]) return(2)                                                                                    # Hetrozygous *H*
      if(x[1] == x[2] && x[1] == snpalleles[1]) return(1)                                                           # Homozygous Allele 1 *A*
      if(x[1] == x[2] && x[1] == snpalleles[2]) return(3)                                                           # Homozygous Allele 2 *B*
      stop("Should not get here, contact Danny.Arends@gmail.com")
    })
    genotypes[snp,] <- genotype
  }

  outCSVR <- rbind(c("Pheno", "", "", peddata[,"Pheno"]), cbind(mapdata[,c("ID","Chr","cM")], genotypes))           # Create CSVRotated output
  write.table(outCSVR, file = out, row.names=FALSE, col.names=FALSE,quote=FALSE, sep=",")                           # Save it to a file
  require(qtl)
  return(read.cross(file=out, "csvr", genotypes=c(1,2,3)))                                                          # Load it using R/qtl read.cross
}

library(qtl)

cross <- convertToCSVR()
cross <- jittermap(cross)
plot(scanone(cross))
