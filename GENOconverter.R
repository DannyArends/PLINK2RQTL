# Convert GeneNetwork .geno files to R/qtl csvr format
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

#Set this to where YOU stored the example data
setwd("D:/Github/PLINK2RQTL/test")

### Convert GeneNetworks .geno format to the R/qtl .csvr input format ###
GENOtoCSVR <- function(genotypes = "BXD.geno", out = "cross.csvr", phenotype = NULL, sex = NULL, verbose = FALSE){
  genodata <- read.table("BXD.geno", sep="\t", skip = 6, header=TRUE, na.strings="U", colClasses="character")
  if(is.null(phenotype)) phenotype <- runif((ncol(genodata)-4))                                                     # If there isn't a phenotype, generate a random one
  if(is.null(sex)) sex <- rep("m", (ncol(genodata)-4))                                                              # If there isn't a sex phenotype, treat all as males
  outCSVR <- rbind(c("Pheno", "", "", phenotype),                                                                   # Phenotype
                   c("sex", "", "", sex),                                                                           # Sex phenotype for the mice
                   cbind(genodata[,c("Locus","Chr", "cM")], genodata[, 5:ncol(genodata)]))                          # Genotypes
  write.table(outCSVR, file = out, row.names=FALSE, col.names=FALSE,quote=FALSE, sep=",")                           # Save it to a file
  require(qtl)
  return(read.cross(file=out, "csvr", genotypes=c("B", "H", "D")))                                                  # Load it using R/qtl read.cross  
}

library(qtl)

cross <- GENOtoCSVR()                                                                                               # Test the conversion from GENENETWORK to R/qtl
cross <- jittermap(cross)
plot(scanone(cross))
