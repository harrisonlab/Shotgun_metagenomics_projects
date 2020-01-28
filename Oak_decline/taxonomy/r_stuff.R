library(tidyverse)
library(data.table)
dat <- fread("LANGDALE.names.out",fill=T,sep="\t")
dat[,c("bin","contig"):=tstrsplit(V2,".fa.",fixed=T)]
dat[,acc:=sub(",.*","",V6)]
dat[,c("kingdom","phylum","class","order","family","genus","species"):=tstrsplit(V8,";",fixed=T)]

prot <- fread("LANGDALE.prots.out",header=F)
setnames(prot,c("acc","protein"))
prot <- unique(prot)
dat <- prot[dat,on=c("acc==acc")]

dat[,(5:10):=NULL]

setnames(dat,c("V1","V2"),c("assigned","fullname"))

dat[,.N,by=c("bin")]
dat[,.N,by=c("bin","kingdom")]

test <- dat[bin=="bin.100",]
test[,.N,by=c("phylum","order")]

countData <- fread("sub_bin.countData")
