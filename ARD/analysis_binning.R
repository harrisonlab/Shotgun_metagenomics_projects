#===============================================================================
#       Load libraries
#===============================================================================
library(DESeq2)
library(BiocParallel)
library(data.table)
library(tidyverse)
# library(devtools)
# install_github("eastmallingresearch/Metabarcoding_pipeline/scripts") (RUN ONCE)
library(metafuncs)
register(MulticoreParam(2)) # watch this as each process will take a copy the dds object - could use a lot of memory for sub_bins
# MultiCore actually shares the memory between parent and worker processes. But, unfortunatly the memory gets copied on R scheduled garbage collection. Annoying
#register(SnowParam(12))# No better - will just have to run on a high memory node (needs ~50G for Langdale sub_bin processing - with multiprocessing)

SITE <- "COMB"
#===============================================================================
#       Load data
#===============================================================================
### count matrix (bins)
countData <- fread(paste0(SITE,".countData"))
setnames(countData,names(countData),sub("(_ND.*_L)([0-9]*)(.*)","_\\2",names(countData)))
setnames(countData,"Domain","PFAM_NAME")

# get number of reads mapped to domains       
print(as.data.frame(colSums(countData[,-1])))

### sub bins
countData <- fread(paste0(SITE,".sub_bins.countData"))
setnames(countData,names(countData),sub("_L","_",names(countData)))
countData[,PFAM_NAME:=gsub("(k[0-9]+_[0-9]+_)(.*)(_[0-9]+_[0-9]+$)","\\2",SUB_BIN_NAME)]

# read in pfam annotation
annotation <- fread("~/pipelines/common/resources/pfam/names.txt")
setnames(annotation,"NAME","PFAM_NAME")
pfam_go <- fread("~/pipelines/common/resources/mappings/pfam_go_map",header=F)

# map bins to pfam and go terms
mapping_pfam <- annotation[countData[,1,with=F],on="PFAM_NAME"]
mapping_pfam <- annotation[countData[,c(1,ncol(countData)),with=F],on="PFAM_NAME"]

mapping_go   <- copy(mapping_pfam)
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
#       Subset DNA/RNA
#===============================================================================
# ensure colData rows and countData columns have the same order
colData <- colData[Sample==names(countData)[-1],]
# remove low count and control samples
myfilter <- colData$Type=="DNA"
# apply filter
colData2 <- copy(colData[which(myfilter==T),])
countData2 <- copy(countData[,c(T,myfilter),with=F])

row_names <- 1 
# creates a dds object and also subsets and orders colData by colnames of countData
dds <- DESeqDataSetFromMatrix(dt_to_df(countData2,row_names), dt_to_df(colData2)[names(countData2)[-row_names],], ~1)
# get size factors
sizeFactors(dds) <-sizeFactors(estimateSizeFactors(dds))

#===============================================================================
#       Differential analysis
#===============================================================================
      
# filter out low counts (especially for subbins - there will be a lot with low counts)
# dds <- dds[rowSums(counts(dds, normalize=T))>5,]

# p value for FDR cutoff
alpha <- 0.1

# the full model
design <- ~Pair+Status 

# set any columns used in model to be factors (deseq should really do this internally...)
dds$Status <- as.factor(dds$Status)
dds$Pair <- as.factor(dds$Pair)

# add full model to dds object
design(dds) <- design

# calculate fit - parallel only useful for subbins
dds <- DESeq(dds,parallel=T)

# calculate results for default contrast (S vs H) - parallel only useful for subbins
res <- results(dds,alpha=alpha,parallel=T)
# merge results with annotation
res_merge <- data.table(inner_join(data.table(SUB_BIN_NAME=rownames(res),as.data.frame(res)),mapping_pfam))
res_merge <- data.table(inner_join(data.table(PFAM_NAME=rownames(res),as.data.frame(res)),annotation))

fwrite(res_merge,paste0(SITE,"_bins.txt"),sep="\t",quote=F)
fwrite(res_merge[padj<=0.1,],paste0(SITE,"_bins_sig.txt"),sep="\t",quote=F)

#===============================================================================
#       Functional analysis (sub bins) 
#===============================================================================
library(topGO)

res_filt <- data.table(left_join(data.table(SUB_BIN_NAME=rownames(res),as.data.frame(res)),mapping_go))
res_filt <- res_filt[complete.cases(res_filt),]
fwrite(res_filt[,toString(V4),by=list(SUB_BIN_NAME)],"topgo_temp",sep="\t",row.names=F,col.names=F,quote=F)
geneID2GO <- readMappings("topgo_temp")
genes <- unique(res_filt[,c(1,3,7)]) # sub_bin_name,fc,p(adjusted)
geneList <- setNames(genes$padj*sign(genes$log2FoldChange),genes$SUB_BIN_NAME)
geneSel <- function(X)abs(X)<=0.05

GOdata <- new("topGOdata",ontology = "BP",allGenes = geneList,geneSel = geneSel,annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)

x="Over"#x="under"
geneSelectionFun(GOdata) <- function(X)abs(X)<=0.05&X>0 # increased expression
geneSelectionFun(GOdata) <- function(X)abs(X)<=0.05&X<0 # decreased expression

# weighted uses the go topology to infer statistical significance
# ks uses strength of signal (i.e. p vlaue) for calculating GO term significance (I'm not convinced this is a good idea)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
#resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultFisher.weight <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher, fisherWeighted=resultFisher.weight,orderBy = "fisherWeighted", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))
#allRes <- GenTable(GOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim, ksWeighted=resultKS.weight,orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))
# over_expressed <- allRes[((allRes$Significant)/(allRes$Expected))>1,]

fwrite(allRes,paste(SITE,x,"GO_RES.txt",sep="_"),sep="\t",quote=F)

pdf(paste(SITE,x,"_GO_plots.pdf",sep="_"))
  showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
 # showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
 # showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
  showSigOfNodes(GOdata, score(resultFisher.weight), firstSigNodes = 5, useInfo = 'all')
dev.off()
#===============================================================================
#       Plots and etc.
#===============================================================================

mypca <- des_to_pca(dds)
d <-t(data.frame(t(mypca$x)*mypca$percentVar))

# Anova of first 4 PC scores
lapply(seq(1:4),function(x) {summary(aov(mypca$x[,x]~Block_pair+Status,colData(dds)))})

# sum of Sum-of-squares 
sum_squares <- t(apply(mypca$x,2,function(x) 
  t(summary(aov(x~Block_pair+Status,colData(dds)))[[1]][2]))
)
colnames(sum_squares) <- c("Block","Condition","residual")
x<-t(apply(sum_squares,1,prop.table))
perVar <- x * mypca$percentVar
#colSums(perVar)
colSums(perVar)/sum(colSums(perVar))*100

# plot with lines joining blocks/pairs
ggsave(paste0(SITE,"_bins.pdf"),plotOrd(d,colData(dds),design="Status",shape="Block_pair",pointSize=2,alpha=0.75,cbPalette=T) )
