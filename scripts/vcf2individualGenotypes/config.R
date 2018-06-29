## BASE OF PROJECT
PROJECT_FOLDER <- "/home/thibaut/Gallus/Projects/genotypeTIGER/"

## ALL THE SNPs (TODO: Have to be sorted?)
#vcf.file.f2 <- paste0(PROJECT_FOLDER, "data/180208.all.223+700.f2.P60_galgal6.vcf.gz")
#vcf.file.f2 <- paste0(PROJECT_FOLDER, "data/180208.all.223+700.f2.P60.vcf.bgz")
vcf.file.f2 <- "/home/thibaut/Gallus/Projects/genotype-F2/Data/resequenced.merged.vcf.bgz"
## THE PEDIGREE
ped1 <-  read.table(paste0(PROJECT_FOLDER,  "data/Ped.f2.f2.f0.txt"),stringsAsFactors = FALSE, header = TRUE, sep = "\t")
ped2 <-  read.table(paste0(PROJECT_FOLDER, "data/Ped.f2.f2.AIL.f0.txt"),stringsAsFactors = FALSE, header = TRUE, sep = "\t")
#f2ped.file <-  paste0(PROJECT_FOLDER, "data/Ped.f2.f2.f0.txt")
f2.ped <- rbind.data.frame(ped1,ped2)
## FILE CONTAINING THE NAMES OF SCAFFOLDS ASSOCIATED WITH A NUMBER
##matchingNames <- paste0(PROJECT_FOLDER, "data/chr_id.galgal6.match.txt")
matchingNames <- paste0(PROJECT_FOLDER, "data/chr_id.match.txt")
#matchingNames <- paste0(PROJECT_FOLDER, "test.txt")
## OUTPUTS: Folder for the intermediate + basename for the final outputs
outFolder <- paste0(PROJECT_FOLDER, "data/with.fam.f2-resequenced/")
cutoffLevel <- 10
nbIndiv <- dim(f2.ped)[1]
fileHap1 <- paste0(outFolder,  "/F2_", nbIndiv, ".within.fam.bins",  cutoffLevel, "_",  Sys.Date())
CORES = 5
## ANNEX SCRIPTS
source(paste0(PROJECT_FOLDER, "scripts/Functions.AIL.impute.R"))
py.script <- paste0("python3 ", PROJECT_FOLDER, "scripts/read.vcf.grand.p.py ") #/usr/local/Cellar/python3/3.6.3/bin/
## PHENOTYPES
phenoFile <- "/mnt/bbg/AIL_reseq/useful_info/phenotypes_f2.csv"
mypheno <-"BW8"
sexFile <- "/mnt/bbg/AIL_reseq/useful_info/AIL_Pedigree_Sex.csv"
mysex <- "Sex"
familyFile <- paste0(PROJECT_FOLDER,  "data/familiesNunber.csv")
myfam <- "fam"
