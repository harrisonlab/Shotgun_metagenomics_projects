library(data.table)
library(tidyverse)
library(metafuncs)

filt <- "grep BPD_transp_2"
files <- list.files(".",".*tab",full.names=F)
qq <- lapply(files,function(x) fread(paste(filt,x)))

names <- sub("\\.tab","",files)
lapply(seq(1:length(qq)),function(i) setnames(qq[[i]],"V2",names[i]))
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qq) 
unsetNA(countData)

countData[,C112:=C112_L1+C112_L2]
countData[,C120:=C120_L1+C120_L2]
colToDelete <- c("C112_L1","C112_L2","C120_L1","C120_L2")
countData[,(colToDelete):=NULL]

apply(countData[,-1],2,function(x) sum(x>0))



# or in shell
FILT=ACR_tran
for FILE in *.tab 
do
  COUNT=$(grep $FILT $FILE|wc -l)
  echo -e "$FILE\t$COUNT"
done



