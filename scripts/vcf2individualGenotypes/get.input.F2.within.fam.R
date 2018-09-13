## Part of Yanjun's script

################################################################################
#############  REQUIRE
################################################################################

suppressWarnings(library(dtplyr, quietly=TRUE))
library(dplyr)
library(argparse, quietly=TRUE)
suppressWarnings(library(pbmcapply, quietly=TRUE)) # Add an overhead on parallel to have a progess bar
##library(parallel)


################################################################################
#############  VARIABLE CREATION
################################################################################


parser <- ArgumentParser()

## specify our desired options
## Needed from config before py.script, vcf.file.f2, f2.ped, matchingNames, outFolder, 

## by default ArgumentParser will add an help option 
parser$add_argument("-p", "--pyscript", type="character", default="read.vcf.grand.p.py", help="Location of the python Script [default \"%(default)\"]")
parser$add_argument("-c", "--cores", type="integer", default=5, help="Number of cores to use [default \"%(default)\"]")
parser$add_argument("pedigree", type="character", help="")
parser$add_argument("vcf", type="character", help="")
parser$add_argument("names", type="character", help="")
parser$add_argument("outFolder", type="character", help="")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

print(args)
pedigree <- args$pedigree
vcf.file.f2 <- args$vcf
outFolder <- args$outFolder
py.script <- args$pyscript
matchingNames <- args$names
CORES <- args$cores

################################################################################
#############  FUNCTION DEFINITION
################################################################################

read_grand.p <- function(pedigreeTable, vcf.file, py, pathout = NULL, generate.input = TRUE, names = "", CORES = 5) {
    cat("1.Please make sure the input F_id_h,M_id_h are from HWS and F_id_l,M_id_l are from LWS","\n")
    ## 1. check if the ped is right
    ## not implemented
    ## 2. read in data
    cat("2.Format data using python","\n")
    if(is.null(pathout)) {
        pathout <- "test.vcf"
    }
    if(generate.input) {
        tmp = group_by(pedigreeTable, fa.h, ma.h, ma.l, fa.l) %>% summarise()
        getFamily <- function(i){
            val <- filter(pedigreeTable, fa.h == tmp[i,]$'fa.h' & ma.h == tmp[i,]$'ma.h' & ma.l == tmp[i,]$'ma.l' & fa.l == tmp[i,]$'fa.l')
            of_id <- val$id.f2
            F_id_h <- paste(".*", tmp[i,]$fa.h, ".*", sep = "")
            M_id_h <- paste(".*", tmp[i,]$ma.h, ".*", sep = "")
            F_id_l <- paste(".*", tmp[i,]$fa.l, ".*", sep = "")
            M_id_l <- paste(".*", tmp[i,]$ma.l, ".*", sep = "")
            return(list("of_id" = of_id, "F_id_h" = F_id_h, "M_id_h" = M_id_h, "F_id_l" = F_id_l, "M_id_l" = M_id_l))
        }
        families = lapply(1:nrow(tmp), getFamily)
        #print(families)
        callScript <- function(family){
            namescmd <- ""
            if (nchar(names) > 0){
                namescmd <- paste0(" --names ", names)
            }
            cmd <- paste(py, namescmd," --ID_of", paste(family$of_id, collapse = " "), "--", family$F_id_h, family$M_id_h, family$F_id_l, family$M_id_l, vcf.file, pathout,sep = " ")
            print(cmd)
            system(cmd)
        }
        pbmclapply(families, callScript, mc.cores = getOption("mc.cores", CORES))
    }
}


################################################################################
#############  EXECUTION
################################################################################

dir.create(outFolder)
f2.ped <-  read.table(pedigree,stringsAsFactors = FALSE, header = TRUE, sep = "\t")
f2.ped$id.f2 = trimws(f2.ped$id.f2)
read_grand.p(f2.ped, vcf.file = vcf.file.f2, pathout = outFolder, py = py.script, generate.input = TRUE, names = matchingNames, CORES=CORES)
