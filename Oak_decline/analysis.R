library(data.table)
library(tidyverse)
library(DESeq2)

qq <- lapply(list.files(".","*",full.names=T),function(x) fread(x))
names<- sub("\\..*","",list.files(".","*",full.names=F,recursive=F))             
qq <- lapply(seq(1:length(qq)),function(i) {X<-qq[[i]];colnames(X)[2] <- names[i];return(X)})

X <- qq %>% reduce(full_join) # plyr method (returns tibble)
X <- Reduce(function(...) merge(..., all = TRUE), qq) # data table method (returns data table)
            
