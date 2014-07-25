# Convert PLINK .ped & .map files to R/qtl csvr format
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

#Set this to where YOU stored the example data
setwd("D:/Github/PLINK2RQTL/test")

### Convert PLINKs .ped and .map files to the R/qtl .csvr input format ###
PLINKtoCSVR <- function(ped = "test_complete.ped", map = "test.map", out = "cross.csvr", missing.genotype = "0", 
                        no.fid = FALSE, no.parents = FALSE, no.sex = FALSE, no.pheno = FALSE,
                        verbose = FALSE){
  mapdata <- read.table(map, colClasses=c("character"))                                                             # Read the MAP data
  colnames(mapdata) <- c("Chr", "ID", "cM", "BP")

  SNPcolnames <- paste0(unlist(lapply(mapdata[,"ID"],rep,2)), c(".A",".B"))                                         # Create column names for the SNPs

  peddata <- scan(ped, what=character(), na.strings=c("-9"))                                                        # Read the PED data, using the much faster scan function
  columnNames <- NULL
  columnNames <- c(columnNames, "IID")
  if(!no.fid)     columnNames <- c(columnNames, "FID")
  if(!no.parents) columnNames <- c(columnNames, "PID", "MID")
  if(!no.sex)     columnNames <- c(columnNames, "Sex")
  if(!no.pheno)   columnNames <- c(columnNames, "Pheno")
  peddata <- matrix(peddata, ncol=length(c(columnNames, SNPcolnames)), byrow = TRUE)

  colnames(peddata) <- c(columnNames, SNPcolnames)

  if(no.pheno) peddata <- cbind(peddata, Pheno=runif(nrow(peddata)))                                                # If there is no phenotype, create a random one
  if(no.sex) peddata <- cbind(peddata, Sex=rep(1, nrow(peddata)))                                                   # If there is no sex, then everyone is a male

  peddata[peddata[,"Sex"] == 1, "Sex"] <- "m"; peddata[peddata[,"Sex"] == 2, "Sex"] <- "f"                          # R/qtl uses m and f, for males and females
  genotypes <- matrix(NA, length(mapdata[,"ID"]), nrow(peddata))                                                    # Empty genotype matrix
  rownames(genotypes) <- mapdata[,"ID"]
  column <- length(columnNames)+1
  for(snp in mapdata[,"ID"]){
    cols <- column:(column+1)                                                                                       # Get the columns associated with this SNP
    snpalleles <- sort(unique(unlist(as.character(peddata[,cols]))))                                                # The SNP alleles
    if(missing.genotype %in% snpalleles){
      snpalleles <- snpalleles[-which(snpalleles == missing.genotype)]                                              # Missing data should not count as an allele
    }
    if(length(snpalleles) > 2){
      cat("[WARNING]", snp, "found multi allelic marker:", snpalleles, ", passed as all missing\n")
      genotype <- rep(NA, nrow(peddata))
    }else{
      if(verbose){
        cat((column - (length(columnNames)+1)) / 2,"/", length(mapdata[,"ID"]), snp,"found", snpalleles,"\n")       # Some info / debug
      }
      genotype <- apply(peddata[,cols], 1, function(x){
        if(x[1] == missing.genotype) return(NA)                                                                     # Missing genotype data
        if(x[1] != x[2]) return(2)                                                                                  # Heterozygous *H*
        if(x[1] == x[2] && x[1] == snpalleles[1]) return(1)                                                         # Homozygous Allele 1 *A*
        if(x[1] == x[2] && x[1] == snpalleles[2]) return(3)                                                         # Homozygous Allele 2 *B*
        stop("Should not get here, contact Danny.Arends@gmail.com")
      })
    }
    column <- column + 2
    genotypes[snp,] <- genotype
  }

  outCSVR <- rbind(c("Pheno", "", "", peddata[,"Pheno"]), 
                   c("sex", "", "", peddata[,"Sex"]),                                                               # Add the sex phenotype
                   cbind(mapdata[,c("ID","Chr","cM")], genotypes))                                                  # Create CSVRotated output
  write.table(outCSVR, file = out, row.names=FALSE, col.names=FALSE,quote=FALSE, sep=",")                           # Save it to a file
  require(qtl)
  return(read.cross(file=out, "csvr", genotypes=c(1,2,3)))                                                          # Load it using R/qtl read.cross
}

library(qtl)

cross <- PLINKtoCSVR("test_complete.ped")                                                                           # Test the conversion from PLINK to R/qtl
plot(scanone(cross))

cross <- PLINKtoCSVR("test_no-pheno.ped", no.pheno=TRUE)                                                            # Test the conversion from PLINK to R/qtl, no phenotype
plot(scanone(cross))

cross <- PLINKtoCSVR("test_minimal.ped", no.fid = TRUE, no.parents = TRUE, no.sex = TRUE, no.pheno = TRUE)          # Test the conversion from PLINK to R/qtl, only IID column
plot(scanone(cross))
