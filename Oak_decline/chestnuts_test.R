library(data.table)
library(tidyverse)
library(metafuncs)

filt <- "grep DDE_Tnp_1_4"
files <- list.files(".",".*tab",full.names=F)
qq <- lapply(files,function(x) {fread(paste(filt,x))})

names <- sub("\\.tab","",files)
invisible(lapply(seq(1:length(qq)),function(i) {setnames(qq[[i]],"V2",names[i])}))
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qq) 

unsetNA(countData)

countData[,C112:=C112_L1+C112_L2]
countData[,C120:=C120_L1+C120_L2]
colToDelete <- c("C112_L1","C112_L2","C120_L1","C120_L2")
countData[,(colToDelete):=NULL]
setnames(countData,sub("_L[0-9]","",names(countData)))

# no. of none zero bins
apply(countData[,-1],2,function(x) {sum(x>0)})

# plot cummulative reads
DT_S <- data.table(Symptom=rowMeans(countData[,as.data.frame(colData[V6=="Symptom","V2",with=F])[,1],with=F]))
DT_H <- data.table(Healthy=rowMeans(countData[,as.data.frame(colData[V6=="Healthy","V2",with=F])[,1],with=F]))
DT <- melt(data.table(DT_H,DT_S))
DT[,value:=log10(value)]
ggplot(DT, aes(x=value, colour = variable)) + stat_ecdf()


# plot something else
DT <- data.table(
          Symptom=rowMeans(countData[,as.data.frame(colData[V6=="Symptom","V2",with=F])[,1],with=F]),
          Healthy=rowMeans(countData[,as.data.frame(colData[V6=="Healthy","V2",with=F])[,1],with=F])
)
DT <- melt(DT)

cut_off=10
ggplot(DT[value>=cut_off,],aes(x=value, fill=variable)) + geom_density(alpha=0.25) + theme_blank()
ggplot(DT[value>=cut_off,],aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
ggplot(DT[value>=cut_off,],aes(x=variable, y=value, fill=variable)) + geom_boxplot()
          
      
colData <- fread("grep Chestnuts colData")
colData <- colData[!duplicated(V2),]

      
# or in shell
FILT=ACR_tran
for FILE in *.tab 
do
  COUNT=$(grep $FILT $FILE|wc -l)
  echo -e "$FILE\t$COUNT"
done
