#===============================================================================
#       Load libraries
#===============================================================================
library(DESeq2)
library(BiocParallel)
library(data.table)
library(tidyverse)
library(devtools)
load_all("~/pipelines/metabarcoding/scripts/myfunctions") # install_github(https://.../)
register(MulticoreParam(12)) # watch this as each process will take a copy the dds object - could use a lot of memory for sub_bins

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

### count matrix (bins)
countData <- fread("countData")
setnames(countData,names(countData),sub("(_ND.*_L)([0-9]*)(.*)","_\\2",names(countData)))
setnames(countData,"Domain","PFAM_NAME")


### sub bins
countData <- fread("countData.subbins")
setnames(countData,names(countData),sub("_L","_",names(countData)))
countData[,PFAM_NAME:=gsub("(k[0-9]+_[0-9]+_)(.*)(_[0-9]+_[0-9]+$)","\\2",SUB_BIN_NAME)]

# read in pfam annotation
annotation <- fread("~/pipelines/common/resources/pfam/names.txt")
setnames(annotation,"NAME","PFAM_NAME")
pfam_go <- fread("~/pipelines/common/resources/mappings/pfam_go_map",header=F)

# drop exact duplicates (subbins only) - not needed
#countData$NAME <- sub("_clust0.*","",countData$V1)
#countdata <- countData[!duplicated(countData[,-1]),-"NAME"]
# countData <- countData[!duplicated(countData[,-1]),]
# countData$BIN_ID <- paste0("BIN",seq(1,nrow(countData)))

# map bins to pfam and go terms
mapping_pfam <- annotation[countData[,c(1,ncol(countData)),with=F],on="PFAM_NAME"]

mapping_go   <- mapping_pfam
mapping_go$ACC <- sub("\\..*","",mapping_go$ACC)
mapping_go <- data.table(left_join(mapping_go,pfam_go,by=(c("ACC"="V1"))))
mapping_go <- mapping_go[complete.cases(mapping_go),]
#bingo_out <- data.table(V1="EMPTY",V2=mapping_go[,BIN_ID],V3=mapping_go[,BIN_ID],V4=NA,V5=mapping_go[,V4],V6=NA,V7="ISS",V8="UNKNOWN",V9="C",V10="UNKNOWN",V11=mapping_go[,BIN_ID],V12="gene",V13="taxon:000",V14="20180101",V15="GD")
#write.table(bingo_out,"gene_association.GO_XXX",col.names=F,row.names=F,quote=F,na="",sep="\t")

# set NA values to 0
#countData[is.na] <- 0
unsetNA(countData)

# remove unused columns from countData (BIN_ID can be mapped to mapping_pfam)
#countData <- countData[,-c("V1","NAME"),with=F]
# sub_bins only
colsToDelete <- c("V1","NAME","PFAM_NAME") 
countData[, (colsToDelete) := NULL]

# read in metadata
colData   <- fread("colData")

#===============================================================================
#       Pool Data/subsample
#===============================================================================
# row_names column of countData object
row_names <- 1 

# creates a dds object and also subsets and orders colData by colnames of countData
dds <- DESeqDataSetFromMatrix(dt_to_df(countData,row_names), dt_to_df(colData)[names(countData)[-row_names],], ~1)

# get number of samples per tree
#sample_numbers <- table(sub("[A-Z]$","",dds$Sample))

# collapse (mean) samples - could just use sum, then sizeFactors will correct
dds <- collapseReplicates(dds,groupby=dds$Sample)

# get size factors
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Differential analysis
#===============================================================================

# remove unhelpful counts (Langdale)
dds <- dds[,dds$Status!="Sandra"]

# filter out low counts (especially for subbins - there will be a lot with low counts)
# dds <- dds[rowSums(counts(dds, normalize=T))>5,]

# p value for FDR cutoff
alpha <- 0.1

# the full model
design <- ~Block_pair+Status # or
design <- ~Status

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
res_merge <- data.table(inner_join(data.table(SUB_BIN_NAME=rownames(res),as.data.frame(res)),mapping_pfam))
res_merge <- data.table(inner_join(data.table(PFAM_NAME=rownames(res),as.data.frame(res)),annotation))

fwrite(res_merge,"CHESTNUTS_pairs.txt",sep="\t",quote=F)
fwrite(res_merge[padj<=0.1,],"CHESTNUTS_pairs_sig.txt",sep="\t",quote=F)

# Langdale
design <- ~Block_pair+Status # no results for ~Status model 
LHS <- c("AOD","COD","Remission","AOD","COD","AOD")
RHS <- c("Healthy","Healthy","Healthy","Remission","Remission","COD")

res <- lapply(seq_along(LHS),function(i) {
  results(dds,alpha=alpha,parallel=T,contrast=c("Status",LHS[i],RHS[i]))
})

invisible(lapply(res,summary))

res_merge <- lapply(res,function(res) {
  data.table(inner_join(data.table(PFAM_NAME=rownames(res),as.data.frame(res)),annotation))
})

lapply(seq_along(res_merge),function(i) {
  invisible(fwrite(res_merge[[i]],paste("LANGDALE",LHS[i],RHS[i],"txt",sep="."),sep="\t",quote=F))
})

#===============================================================================
#       Functional analysis
#===============================================================================
library(topGO)
res_filt <- data.table(left_join(data.table(BIN_ID=rownames(res),as.data.frame(res)),mapping_go))
res_filt <- res_filt[complete.cases(res_filt),]
write.table(res_filt[,toString(V4),by=list(BIN_ID)],"topgo_temp",sep="\t",row.names=F,col.names=F,quote=F)
geneID2GO <- readMappings("topgo_temp")
genes <- unique(res_filt[,c(1,3,7)])
geneList <- setNames(genes$padj*sign(genes$log2FoldChange),genes$BIN_ID)
geneSel <- function(X)abs(X)<=0.05

GOdata <- new("topGOdata",ontology = "BP",allGenes = geneList,geneSel = geneSel,annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)

geneSelectionFun(GOdata) <- function(X)abs(X)<=0.05&X>0 # increased expression
geneSelectionFun(GOdata) <- function(X)abs(X)<=0.05&X<0 # decreased expression


resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

resultKS.weight <- runTest(GOdata, algorithm = "weight01", statistic = "ks")

allRes <- GenTable(GOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))
over_expressed <- allRes[((allRes$Significant)/(allRes$Expected))>1,]


pdf("plot.pdf")
  showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
dev.off()
#===============================================================================
#       Plots and etc.
#===============================================================================

mypca <- des_to_pca(dds)
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

plotOrd(d,colData(dds),design="Status")
