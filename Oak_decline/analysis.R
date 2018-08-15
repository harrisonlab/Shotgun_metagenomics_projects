library(data.table)
library(dtplyr)

qq <- lapply(list.files(".","*",full.names=T),function(x) fread(x))
             
X <-              
