#===============================================================================
#       Load libraries
#===============================================================================
library(DESeq2)
library(BiocParallel)
library(data.table)
library(tidyverse)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions") # install_github(https://.../)
register(MulticoreParam(12))

#===============================================================================
#       Load data
#===============================================================================

####  counts tables
# get list of count tables
qq    <- lapply(list.files(".","*",full.names=T),function(x) {fread(x)})

# get count file names and substitute to required format
names <- sub("(_ND.*_L)([0-9]*)(.*)","_\\2",list.files(".","*",full.names=F,recursive=F))

# apply names to appropriate list columns (enables easy joining of all count tables)
qq    <- lapply(seq(1:length(qq)),function(i) {X<-qq[[i]];colnames(X)[2] <- names[i];return(X)})

# merge count tables into single table
countData <- Reduce(function(...) {merge(..., all = TRUE)}, qq) # data table method (returns data table)
#countData    <- qq %>% purr::reduce(full_join,by="V1") # plyr method (returns tibble - but not always...)

### count matrix
countData <- fread("countData")
names(countData) <- sub("(_ND.*_L)([0-9]*)(.*)","_\\2",names(countData))

### sub bins
countData <- fread("countData.subbins")

# drop exact duplicates (subbins only)
countData$pfam <- sub("_clust0.*","",countData$V1)
countdata <- countData[!duplicated(countData[,-1]),-"pfam"]

# set NA values to 0
countData[is.na] <- 0

# read in metadata
colData   <- fread("colData")

# subset metadata 
colData <- colData[SampleID%in%names(countData),]

# read in pfam annotation
annotation <- fread("~/pipelines/common/resources/pfam/names.txt")

#===============================================================================
#       Pool Data/subsample
#===============================================================================

dds <- DESeqDataSetFromMatrix(dt_to_df(countData), dt_to_df(colData), ~1)

# get number of samples per tree
#sample_numbers <- table(sub("[A-Z]$","",dds$Sample))

# collapse (mean) samples - could just use sum, then sizeFactors will correct
dds <- collapseReplicates2(dds,groupby=dds$Sample,simple=T)

# set the dds sizefactor to the number of samples
# dds$sizeFactor <- as.vector(sample_numbers/2)

# recreate countData and colData
#countData<- round(counts(dds,normalize=T),0)
#colData <- as.data.frame(colData(dds))

# new dds object with the corrected data set
#dds <- DESeqDataSetFromMatrix(countData,colData,~1)

#===============================================================================
#       Differential analysis
#===============================================================================

# create dds object
# dds <- DESeqDataSetFromMatrix(dt_to_df(countData), dt_to_df(colData), ~1)

# get size factors
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

# filter out low counts (especially for subbins - there will be a lot with low counts)
# dds <- dds[rowSums(counts(dds, normalize=T))>5,]

# p value for FDR cutoff
alpha <- 0.1

# the full model
design <- ~Block_pair+Status

# set any columns used in model to be factors (deseq should really do this internally...)
dds$Status <- as.factor(dds$Status)
dds$Block_pair <- as.factor(dds$Block_pair)

# add full model to dds object
design(dds) <- design

# calculate fit - parallel only useful for subbins
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H) - parallel only useful for subbins
res <- results(dds,alpha=alpha,parallel=T)

# merge results with annotation
res_merge <- data.table(inner_join(data.table(NAME=rownames(res),as.data.frame(res)),annotation))


#===============================================================================
#       Plots and etc.
#===============================================================================

mypca <- des_to_pca(dds)
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

plotOrd(d,colData(dds),design="Status")
