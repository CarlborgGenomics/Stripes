## Part of Yanjun's script

################################################################################
#############  REQUIRE
################################################################################

library(dtplyr)
library(dplyr)

require(pbmcapply) # Add an overhead on parallel to have a progess bar
##require(parallel)

################################################################################
#############  VARIABLE CREATION
################################################################################

source("config.R")

################################################################################
#############  FUNCTION DEFINITION
################################################################################

#source("Functions.AIL.impute.R")
#source("get.input.F2.functions.R")

################################################################################
#############  EXECUTION
################################################################################

dir.create(outFolder)
name <-  "test"
fileHap1 <- paste0(outFolder,  "/F2_", name, "within.fam.bins", Sys.Date())
############# START THE GENERATION OF FIXED VCF ###########
header.f2 <- vcf.header2(vcf = vcf.file.f2) # check the vcf to getout the sample ID. there might be a broken pipe error,  but ignore it
##f2 <- header.f2[grep(pattern = "F2.*", x = header.f2)] ## will be all F2 in the end,  for now a subset since all not available, NOT REALLY GENERIC
f2 <- header.f2[grep(pattern = ".*F2.*", x = header.f2)] 
id2 <- gsub(replacement = "\\1", pattern = "F2_(.*)_(S|m).*", x = f2)
f2.ped$id.f2 = trimws(f2.ped$id.f2)
read_grand.p(f2.ped, vcf.file = vcf.file.f2, pathout = outFolder, py = py.script, generate.input = TRUE, names = matchingNames)
