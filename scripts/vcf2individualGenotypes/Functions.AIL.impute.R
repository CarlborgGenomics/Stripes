
#######################################Functions need in AIL imputed #################################################
## Be aware that some of the functions not menioned in README.R are rather preliminary.

#############>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Read in data <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Argument in this fuction are offspring id, father id, mother id, vcf file name and path for output files, whether generate input or not
read_trio <- function(of_id,F_id,M_id,vcf.file,pathout=NULL,generate.input=T){
  if(!require(data.table))
    require(data.table)
  cat("1.Validte pedigree","\n")
  #1. check if the ped is right
  # not implemented 
  #2. read in data
  cat("2.Format data using python","\n")
  if(is.null(pathout)){
    if(generate.input){
      bash  <- paste(py.script,of_id," ",F_id," ",M_id,vcf.file ," ~/Documents/impute/data/testplate.F1/test.vcf",sep="")
      system(bash)
    }
    cat("3. Read in data using fread")
    input <- fread("~/Documents/impute/data/testplate.F1/test.vcf")
  }else{
    if(generate.input){
    bash  <- paste(py.script,of_id," ",F_id," ",M_id,vcf.file ," ",pathout,sep="")
    system(bash)
    }
    cat("3. Read in data using fread")
    #cat("FAIL HERE")
    #files <- gsub(pattern = "\\s+(.*)", replacement="\\1",pathout)
    #cat("#########\n")
    #cat(files)
    #cat("#########\n")
    input <- fread(gsub(pattern = "\\s+(.*)", replacement="\\1",pathout))
  }
  colnames(input) <- c("chr","pos","ref","alt","f1","f2","m1","m2","of1","of2")
  return(input)
}

############
get.info <- function(chroms,chroms.len,bin.size=1e6){
  num.bin <- rep(NA,length(chroms))
  index <- data.frame(array(NA,dim = c(length(chroms),2)))
  colnames(index) <- c("start","end")
  rownames(index) <- chroms
  num.bin <- ceiling(chroms.len/bin.size)
  names(num.bin) <- chroms
  
  if(any(num.bin==1))
    warnings("a few chr only have one bin, start and end are set to the same","\n")
  
  index$start[1] <- 1
  index$end[1] <- num.bin[1]
  loca <- seq(from = 0.5,to =num.bin[1])
  loca.chr <- rep(chroms[1] , num.bin[1])
  for( i in 2:length(chroms)){
    index$start[i] <- sum(num.bin[1:c(i-1)])+1
    index$end[i] <- sum(num.bin[1:i])
    loca <- c(loca,seq(from = 0.5,to =num.bin[i]))
    loca.chr <- c(loca.chr,rep(chroms[i],num.bin[i]))
  }
  chr.loca <- rep(NA,max(index$end))
  for( i in 1:length(loca)){
    chr.loca[i] <- paste(loca.chr[i],"-",loca[i],sep = "")  
  }
  return(list("num.bin" = num.bin,"index"=index,"loca"=loca,"loca.chr"=loca.chr,"chr.loca"=chr.loca))
}


# # input <-out
# _ped <- function(input){
#   hetsite <- input$of1 !=input$of2 & (input$f1!=input$f2 | input$m1!=input$m2)
#   hetsite <-  (input$f1!=input$f2 & input$m1!=input$m2)
#   
#   info.het <-input[hetsite,]
#   
#   #homo alt par
#   home.alt <- input$f1==input$f2 & input$m1==input$m2
#   info.home.alt <- input[home.alt,]
#   #input[input$f1==input$f2 & input$m1==input$m2 & input$of1 !=input$of2,]
# }
# 
# # input <- out
# Hom_alt <- function(input){
#   hom.alt.site <- (input$f1==input$f2 & input$m1==input$m2 & input$f1 != input$m1)
# 
#   info.het <-input[hom.alt.site,]
#   
#   #homo alt par
#   home.alt <- input$f1==input$f2 & input$m1==input$m2
#   info.home.alt <- input[home.alt,]
#   #input[input$f1==input$f2 & input$m1==input$m2 & input$of1 !=input$of2,]
# }

# input <- read_trio(of_id = "1343_F2_S12_L001",F_id = F_id,M_id = M_id,vcf.file = vcf.file.f2)

########## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> tracing
# input <- data1
tracing_physical <- function(input,bin=1e6,chr.len=200*1e6,cut=5){
    if(!require(zoo))
        require(zoo)
    # We assume f1 == f2 AND f1 != m1 AND m1 == m2
    trace.of1 <- input[,"f1"] - input$of1
    trace.of2 <- input[,"f1"] - input$of2
    trace.of <- abs(trace.of1 + trace.of2) / 2
    num.bin <- chr.len/bin
    bin.loca <- seq(0,chr.len,by=bin)+1
    index.pos <- findInterval(input$pos,bin.loca)
    # Because of High coverage need to check of1 and of2 (if low coverage replace trace.of1 and comment 2 lines)
    out <- calculate.origin.phy(x = trace.of, index.pos = index.pos,cut=cut) 
    f <- as.numeric(out$ratio.mrk)
    #m <- as.numeric(calculate.origin.phy(x = trace.mat1$m1, index.pos = index.pos,cut=cut)$ratio.mrk)
    fn <- as.numeric(out$num.mrk)
    loca <- (sort(unique(index.pos),decreasing = F)*bin - 0.5*bin)
    ##return(list("f"=f,"m"=m,"fn"=fn,"loca"=loca))
    return(list("f"=f, "fn"=fn, "loca"=loca))
}

calculate.origin.phy <- function(x,index.pos=index.pos,cut=5){
    ratio <- function(x) {
        if(length(x) > cut ){
            return(sum(x)/length(x))
        }else{
            return(NA)
        }
    }
    num.mrk <- tapply(X = x, INDEX = index.pos, FUN = function(x) length(x))
    ratio.mrk <- tapply(X = x, INDEX = index.pos, FUN = ratio )
    return(list("num.mrk"=num.mrk,"ratio.mrk"=ratio.mrk))
}

calculate.orgin <- function(x,index=index){
  return(tapply(X = x, INDEX = index, FUN = function(x) sum(x==0,na.rm = T)))
}

calculate.orgin_step <- function(x,stepsize,step){
  
  return(rollapply(x,width=step,by=stepsize,FUN=function(x) sum(x==0,na.rm = T)))
}
# input <- out
tracing <- function(input,step=20,stepsize=NULL){
  if(!require(zoo))
    require(zoo)
  trace.mat1 <- input[,c("f1","f2","m1","m2")] - input$of1
  
  if(is.null(stepsize)){
    total <- ceiling(nrow(trace.mat1)/step)
    total.floor <- floor(nrow(trace.mat1)/step)
    index <- rep(1:total,each=step)[1:nrow(trace.mat1)]
    #loca <- aggregate(data$pos,by=list(index),FUN = mean)/10e5
    #<<<<<<<<<< pos.start pos.end
    cat("huge gaps within win are not checked for now","\n")
    pos <- Bin_pos.nosliding(input = input,index=index,trace.mat1 = trace.mat1,total = total,total.floor = total.floor)
    return(list("counts"= apply(X = trace.mat1, MARGIN = 2, FUN = calculate.orgin,index=index)/step,"pos"=pos))
    #return(list("count"=apply(X = trace.mat1, MARGIN = 2, FUN = calculate.orgin,index=index),"pos"=pos))
  }else{
    count.step <- apply(X = trace.mat1, MARGIN = 2, FUN = calculate.orgin_step,step=step,stepsize=stepsize)/step
    pos <- Bin_pos_sliding(input = input,count.step = count.step,step = step,stepsize = stepsize)
    return(list("counts"= count.step ,"pos"=pos))
  }
  
}


## function to get bin 

Bin_pos.nosliding <- function(input,index,trace.mat1,total,total.floor){
  # check if the bin is equal
  pos <- data.frame(array(NA,dim = c(total,2)))
  colnames(pos) <- c("start","end")
  pos$start <- input$pos[c(0:c(total-1))*step+1]
  
  if(nrow(trace.mat1) %% step == 0){
    cat(" window size is equal","\n")
    pos$end <- input$pos[c(0:c(total-1))*step+step]
    #pos$end[total] <- input$pos[nrow(trace.mat1)]
  }else{
    pos$end <- input$pos[c(0:c(total-1))*step+step]
    pos$end[total] <- input$pos[nrow(trace.mat1)]
  }
  return(pos)
}

##########################Bin sliding 
Bin_pos_sliding <- function(input,count.step,step,stepsize){
  if(stepsize!=1)
    stop("stepsize !=1 is not implemetended yet","\n")
  len <- nrow(count.step)
  pos <- data.frame(array(NA,dim = c(len,2)))
  colnames(pos) <- c("start","end")
  pos$start <- input$pos[1:len]
  pos$end <- input$pos[0:c(len-1)+step]
  return(pos)
}
##### trc to hap matrix
# trc<- tracing(input = data,step = win,stepsize = 1)

Trace2hap <- function(thres=0.85,count=trc$counts,pos = trc$pos,win=NULL){
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  #>>>>>>>>>>>>>>>>>> Marking homozogous region >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(!is.null(win))
    stop("sliding windown is not implemented")
  ff <- Set.hap_nosliding(above = which(trc$counts[,"f1"] >=0.85),pos = pos,line = "ff" )
  mm <- Set.hap_nosliding(above = which(trc$counts[,"m1"] >=0.85),pos = pos,line = "mm" )
  out <- rbind.data.frame(ff$hom,mm$hom,ff$het,mm$het)
  return(out)
}
# 
# 
# 
# win <- 70
# step = win
# trc<- tracing(input = data,step = win,stepsize = 1)
# thres <- 0.85
# above <- which(trc$counts[,"f1"] >=0.85)
# pos =trc$pos
# gap <- ceiling(win - win*(thres-0.5)*2)
# 
# #plot hom
# for( i in 1:nrow(hom)){
#   segments(x0 = as.numeric(hom[i,1])/1e6,y0 =1.05,x1 = as.numeric(hom[i,2])/1e6,col = "blue",cex=2 )
# }
# 
# 
# trc$counts[,"f1"][id.het[1]:id.het[2]]

Set.hap_sliding <- function(above,pos,line="ff",thres=thres,win=win){
  dis <-  pos$start[above][2:length(above)] -pos$end[above][1:(length(above)-1)]
  id <- which(dis >0)
  id.het <- c(above[1],above[id],above[length(above)])
  het <- data.frame(array("Unknow",dim = c(length(id.het),3)),stringsAsFactors = F)
  if(length(unique(id.het))==1){
    het <- c(pos$start[id],pos$start[id])
  }else{
    for( i in 1:length(id.het)){
      het[i,1] <- pos$start[id.het[i]]
      het[i,2] <- pos$end[id.het[i]]
    }
  }
  
  hom <- data.frame(array(line,dim = c(c(length(id)+1),3)),stringsAsFactors = F)
  gap <- ceiling(win - win*(thres-0.5)*2)
  warning(" This only works for sliding window, Recombination break point is estimated in a very rough way ")
  if(length(id)==1){
    hom[1,1] <- pos$start[above[1]+gap]
    hom[1,2] <- pos$start[above[id]-gap]
    hom[2,1] <- pos$start[above[id+1]+gap]
    hom[2,2] <- pos$start[above[length(above)]-gap]
  }else{
    for( i in 1:c(length(id)+1)){
      if(i==1){
        hom[i,1] <- pos$start[above[1]+gap]
        hom[i,2] <- pos$start[above[id[i]]-gap]
      }else if (i >1 & i != c(length(id)+1)){
        
        hom[i,1] <- pos$start[above[id[i-1]+1]+gap]
        hom[i,2] <- pos$start[above[id[i]]-gap]
        
      } else if(i == c(length(id)+1)) {
        hom[i,1] <- pos$start[above[id[i-1]+1]+gap]
        hom[i,2] <- pos$start[above[length(above)]-gap]
      }
    }
  }
  return(list("het"=het,"hom"=hom))
  
}


## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< estimate mrk density
density.mrk <- function(input,len=200*10e6){
  return(len/nrow(input))
}
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> check num of markers in the chromsome
# file <- "~/Documents/impute/data/255/171113.Fixed.sites.founder.all.chr.txt"
# match.file <- fread("~/Documents/impute/data/chr_id.match.txt")
Validate.chr<- function(file){
  if(!require(data.table))
    require(data.table)
  file.in <- fread(file)
  chr <- unique(file.in$V1)
  chr.num <- numeric(length(chr))
  for( i in 1:length(chr)){
    chr.num[i] <- sum(file.in$V1==chr[i]) 
  }
  names(chr.num) <- chr
  return(chr.num)
}
# chr.num <- .chr(file)
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  first bin : check gaps
Validate.geno <- function(input){
  len  <- nrow(input)
  dis <- input$pos[2:len]-input$pos[1:c(len-1)]
  return(dis)
}


######## crap
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>> Marking homozogous region >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# hap <- rep("fm",nrow(count))
# ff <- which(trc$counts[,"f1"] >=0.85)
# mm <- which(trc$counts[,"m1"] >=0.85)
# Set.hap(above = which(trc$counts[,"f1"] >=0.85),pos = pos,line = "ff" )
# hap[ff] <- "ff"
# hap[mm] <- "mm"
# #out <- cbind(pos,hap)
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>> in father haplotype 
# dis <-  pos$start[ff][2:length(ff)] -pos$end[ff][1:(length(ff)-1)]
# id <- which(dis >0)
# homf <- data.frame(array("ff",dim = c(c(length(id)+1),3)),stringsAsFactors = F)
# if(length(id)==1){
#   homf[1,1] <- pos$start[ff[1]]
#   homf[1,2] <- pos$end[ff[id]]
#   homf[2,1] <- pos$start[ff[id+1]]
#   homf[2,2] <- pos$end[ff[length(ff)]]
# }else{
# for( i in 1:c(length(id)+1)){
#   if(i==1){
#     homf[i,1] <- pos$start[ff[1]]
#     homf[i,2] <- pos$end[ff[id[i]]]
#   }else if (i >1 & i != c(length(id)+1)){
#     
#     homf[i,1] <- pos$start[ff[id[i-1]+1]]
#     homf[i,2] <- pos$end[ff[id[i]]]
#     
#   } else if(i == c(length(id)+1)) {
#     homf[i,1] <- pos$start[ff[id[i-1]+1]]
#     homf[i,2] <- pos$end[ff[length(ff)]]
#   }
#   
# }
# }
# # >>>>>>>>>>>>>>>>>>>>>>>>>>>> in mother haplotype 
# dis2 <-  pos$start[mm][2:length(mm)] -pos$end[mm][1:(length(mm)-1)]
# id2 <- which(dis2 >0)
# homm <- data.frame(array("mm",dim = c(c(length(id2)+1),3)),stringsAsFactors = F)
# if(length(id2)==1){
#   homm[1,1] <- pos$start[mm[1]]
#   homm[1,2] <- pos$end[mm[id2]]
#   homm[2,1] <- pos$start[mm[id2+1]]
#   homm[2,2] <- pos$end[mm[length(mm)]]
# }else{
#   for( i in 1:c(length(id2)+1)){
#     if(i==1){
#       homm[i,1] <- pos$start[mm[1]]
#       homm[i,2] <- pos$end[mm[id2[i]]]
#     }else if (i >1 & i != c(length(id2)+1)){
#       
#       homm[i,1] <- pos$start[mm[id[i-1]+1]]
#       homm[i,2] <- pos$end[mm[id[i]]]
#       
#     } else if(i == c(length(id)+1)) {
#       homm[i,1] <- pos$start[mm[id[i-1]+1]]
#       homm[i,2] <- pos$end[mm[length(mm)]]
#     }
#     
#   }
# }
# out <- rbind(homf,homm)  
# colnames(out) <- c("start","end" ,"hap")
# cat("gaps within window is not checked")
# return(out)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>> Marking heterozogous region >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Set.hap_nosliding <- function(above,pos,line="ff"){
  dis <-  pos$start[above][2:length(above)] -pos$end[above][1:(length(above)-1)]
  id <- which(dis >0)
  id.het <- c(above[1],above[id],above[length(above)])
  het <- data.frame(array("Unknow",dim = c(length(id.het),3)),stringsAsFactors = F)
  if(length(unique(id.het))==1){
    het <- c(pos$start[id],pos$start[id])
  }else{
    for( i in 1:length(id.het)){
      het[i,1] <- pos$start[id.het[i]]
      het[i,2] <- pos$end[id.het[i]]
    }
  }
  
  hom <- data.frame(array(line,dim = c(c(length(id)+1),3)),stringsAsFactors = F)
  warning(" This only works for non-sliding window, Recombination break point is included in the borader of hom region, remember to overwrite using het region in python")
  if(length(id)==1){
    hom[1,1] <- pos$start[above[1]]
    hom[1,2] <- pos$end[above[id]]
    hom[2,1] <- pos$start[above[id+1]]
    hom[2,2] <- pos$end[above[length(above)]]
  }else{
    for( i in 1:c(length(id)+1)){
      if(i==1){
        hom[i,1] <- pos$start[above[1]]
        hom[i,2] <- pos$end[above[id[i]]]
      }else if (i >1 & i != c(length(id)+1)){
        
        hom[i,1] <- pos$start[above[id[i-1]+1]]
        hom[i,2] <- pos$end[above[id[i]]]
        
      } else if(i == c(length(id)+1)) {
        hom[i,1] <- pos$start[above[id[i-1]+1]]
        hom[i,2] <- pos$end[above[length(above)]]
      }
    }
  }
  return(list("het"=het,"hom"=hom))
  
}

