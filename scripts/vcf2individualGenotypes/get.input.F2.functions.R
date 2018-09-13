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
    #cat("3. Read in data using fread")
    #readFile <- function(pathFile){
    #    if (file.size(pathFile) > 0){
    #        #print(pathFile)
    #        input <- fread(pathFile)
    #        #colnames(input) <- c("chr","pos","ref","alt","f1","f2","m1","m2","of1","of2")
    #        return(input)
    #    }else {return(NA)}
    #}
    #createPath <- function(x){return(paste0(pathout, "/", x, ".vcf"))}
    #pathList <- lapply(pedigreeTable$id.f2, createPath)
    #pathList <- list.files(pathout, pattern = ".vcf$", full.names = TRUE)
    #inputList <- lapply(pathList, readFile)
    #return(inputList)
}



vcf.header2 <- function(vcf) {
    vcf <- gsub(pattern = "\\s+(.*)",replacement = "\\1",vcf)
    # readLines(gzcon(file("your_file.txt.gz", "rb")))
    con = file(vcf, "r")
    line = readLines(con, n = 1)
    c <- 0
    while ( substring(line[1], 1, 1) == "#" ) {
        c <- c + 1
        line = readLines(con, n = 1)
        if ( length(line) == 0 ) {
            break
        }
        if (any(grep("#CHR", line))){
            b.s <- strsplit(line,"\t")
            break
        }
    }
    close(con)
    #print(c)
    if(class(b.s) == "list") {
        return(unlist(b.s))
    } else {
        return(b.s)
    }
}
