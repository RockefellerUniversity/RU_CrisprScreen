params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE
library(DESeq2)
library(GenomicAlignments)
library(Rsubread)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Crispr Screenining analysis in R

---
"    
  )
  
}



## ----downloadFile,message=FALSE-----------------------------------------------
download.file("https://github.com/LewisLabUCSD/PinAPL-Py/archive/master.zip","pieappleData.zip")
unzip("pieappleData.zip")
download.file("http://pinapl-py.ucsd.edu/example-data","TestData.zip")
unzip("TestData.zip")
download.file("http://pinapl-py.ucsd.edu/run/download/example-run","Results.zip")
unzip("Results.zip")


## ----shortreada,include=FALSE-------------------------------------------------
library(ShortRead)



## ----shortread, warning=F, message=F------------------------------------------
library(ShortRead)


## ----crsiprRep1Reads,echo=T,eval=T--------------------------------------------
fqSample <- FastqSampler("PinAPL-Py-master/Data/Tox-A_R01_98_S2_L008_R1_001_x.fastq.gz",n=10^6)
fastq <- yield(fqSample)


## ----mycRep1ReadsShortReadQ---------------------------------------------------
fastq


## ----mycRep1ReadsAccessor-----------------------------------------------------
readSequences <- sread(fastq)
readQuality <- quality(fastq)
readIDs <- id(fastq)
readSequences


## ----mycRep1ReadsQScores------------------------------------------------------
readQuality <- quality(fastq)
readQualities <- alphabetScore(readQuality)
readQualities[1:10]


## ----mycRep1ReadsQScoresPlot,fig.height=3,fig.width=8-------------------------
library(ggplot2)
toPlot <- data.frame(ReadQ=readQualities)
ggplot(toPlot,aes(x=ReadQ))+geom_histogram()+theme_minimal()


## ----mycRep1ReadsAlpFreq------------------------------------------------------
readSequences <- sread(fastq)
readSequences_AlpFreq <- alphabetFrequency(readSequences)
readSequences_AlpFreq[1:3,]


## ----mycRep1ReadsAlpFreqSum---------------------------------------------------
summed__AlpFreq  <- colSums(readSequences_AlpFreq)
summed__AlpFreq[c("A","C","G","T","N")]


## ----mycRep1ReadsAlpByCycle---------------------------------------------------
readSequences_AlpbyCycle <- alphabetByCycle(readSequences)
readSequences_AlpbyCycle[1:4,1:10]


## ----mycRep1ReadsAlpByCyclePlot-----------------------------------------------
AFreq <- readSequences_AlpbyCycle["A",]
CFreq <- readSequences_AlpbyCycle["C",]
GFreq <- readSequences_AlpbyCycle["G",]
TFreq <- readSequences_AlpbyCycle["T",]
toPlot <- data.frame(Count=c(AFreq,CFreq,GFreq,TFreq),
                     Cycle=rep(1:max(width(readSequences)),4),
                     Base=rep(c("A","C","G","T"),each=max(width(readSequences))))



## ----mycRep1ReadsAlpByCyclePlot2,fig.height=4,fig.width=8---------------------

ggplot(toPlot,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
  theme_bw()


## ----mycRep1ReadsQByCycle-----------------------------------------------------
qualAsMatrix <- as(readQuality,"matrix")
qualAsMatrix[1:2,]


## ----mycRep1ReadsQByCyclePlotfig.width=8,fig.height=4-------------------------
toPlot <- colMeans(qualAsMatrix)
plot(toPlot)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Aligning data

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Aligning data

---
"    
  )
  
}



## ----fa2----------------------------------------------------------------------
GeCKO <- read.delim("PinAPL-py_demo_data/GeCKOv21_Human.tsv")
GeCKO[1:2,]


## ----fa3----------------------------------------------------------------------
require(Biostrings)
sgRNAs <- DNAStringSet(GeCKO$seq)
names(sgRNAs) <- GeCKO$UID


## ----fa4----------------------------------------------------------------------
writeXStringSet(sgRNAs,
                file="GeCKO.fa")


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
library(Rsubread)
buildindex("GeCKO","GeCKO.fa", 
           indexSplit=FALSE)



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs <- "PinAPL-py_demo_data/Control_R1_S14_L008_R1_001_x.fastq.gz"
myMapped <- align("GeCKO",myFQs,output_file = gsub(".fastq.gz",".bam",myFQs),
                  nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA")



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped 



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs <- "PinAPL-py_demo_data/Control_R1_S14_L008_R1_001_x.fastq.gz"
myMapped <- align("GeCKO",myFQs,output_file = gsub(".fastq.gz",".bam",myFQs),
                  nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1)



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped 



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs <- "PinAPL-py_demo_data/Control_R1_S14_L008_R1_001_x.fastq.gz"
myMapped <- align("GeCKO",myFQs,output_file = gsub(".fastq.gz",".bam",myFQs),
                  nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1,
                  maxMismatches = 0,indels = 0)



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped 



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
require(GenomicAlignments)
temp <- readGAlignments(gsub(".fastq.gz",".bam",myFQs))
temp


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
temp


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigars <- cigar(temp)
cigars[1:5]


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigarRLE <- cigarToRleList(cigars)
cigarRLE[1]


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigarMat <- matrix(as.vector(unlist(cigarRLE)),ncol=50,byrow = TRUE)
cigarMat[1:2,]


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigarFreq <- apply(cigarMat,2,table)
cigarFreq[1:2,]


## ---- echo=TRUE,eval=TRUE,fig.height=3,fig.width=8----------------------------
require(ggplot2)
toPlot <- data.frame(Freq=c(cigarFreq["S",],cigarFreq["M",]),
                      Cigar=rep(c("S","M"),each=ncol(cigarFreq)),
                      Cycle=rep(seq(1,ncol(cigarFreq)),2))
ggplot(toPlot,aes(x=Cycle,y=Freq,colour=Cigar))+geom_line()+theme_bw()


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
counts <- data.frame(table(seqnames(temp)),row.names = "Var1")
counts[1:2,,drop=FALSE]


## ---- echo=TRUE,eval=TRUE,fig.height=3,fig.width=8----------------------------
ss <- read.delim("New_Example_1580342908_4/Analysis/01_Alignment_Results/ReadCounts_per_sgRNA/Control_1_GuideCounts.txt")
new <- merge(ss,counts,by.x=1,by.y=0)
plot(new$X0,new$Freq)


## ---- echo=TRUE,eval=TRUE,fig.height=3,fig.width=8----------------------------
counts <- data.frame(table(seqnames(temp[width(temp) == 20])),row.names = "Var1")
new <- merge(ss,counts,by.x=1,by.y=0)
plot(new$X0,new$Freq)


## ---- echo=TRUE---------------------------------------------------------------
load("data/sgRNACounts.RData")
sgRNAcounts[1:4,]


## ---- echo=TRUE---------------------------------------------------------------
library(DESeq2)
metadata <- DataFrame(Group=factor(c("Control","Control","ToxB","ToxB"),
                                   levels = c("Control","ToxB")),
                      row.names = colnames(sgRNAcounts))
dds <- DESeqDataSetFromMatrix(sgRNAcounts,colData = metadata,design = ~Group)


## ---- echo=TRUE---------------------------------------------------------------
dds <- DESeq(dds)


## ---- echo=TRUE,fig.height=3,fig.width=8--------------------------------------
normCounts <- counts(dds,normalized=TRUE)
boxplot(log2(normCounts+0.1))


## ---- echo=TRUE---------------------------------------------------------------
ToxBvsControl <- results(dds,contrast=c("Group","ToxB","Control"))
ToxBvsControl <- ToxBvsControl[order(ToxBvsControl$pvalue),]
ToxBvsControl


## ---- echo=TRUE,fig.height=3,fig.width=8--------------------------------------
ss <- read.delim("New_Example_1580342908_4/Analysis/02_sgRNA-Ranking_Results/sgRNA_Rankings/ToxB_avg_0.01_Sidak_sgRNAList.txt")
toPlt <- merge(ss,as.data.frame(ToxBvsControl),by.x=1,by.y=0)
corr <- cor(log2(toPlt[!is.na(toPlt$padj),]$fold.change),toPlt[!is.na(toPlt$padj),]$log2FoldChange)
plot(log2(toPlt[!is.na(toPlt$padj),]$fold.change),toPlt[!is.na(toPlt$padj),]$log2FoldChange,main=corr)


## -----------------------------------------------------------------------------
ToxBvsControl <- as.data.frame(ToxBvsControl)[order(ToxBvsControl$pvalue),]
ToxBvsControl$Enriched <- !is.na(ToxBvsControl$padj) & ToxBvsControl$pvalue < 0.05 &ToxBvsControl$log2FoldChange > 0
ToxBvsControl[1:2,]


## -----------------------------------------------------------------------------
ToxBvsControl <- merge(GeCKO,ToxBvsControl,by.x=2,by.y=0)
ToxBvsControl[1:2,]


## ----eval=FALSE---------------------------------------------------------------
## genes <- unique(ToxBvsControl$gene_id)
## listofGene <- list()
## for(i in 1:length(genes)){
##   tempRes <- ToxBvsControl[ToxBvsControl$gene_id %in% genes[i],]
##   meanLogFC <- mean(tempRes$log2FoldChange,na.rm=TRUE)
##   logFCs <- paste0(tempRes$log2FoldChange,collapse=";")
##   minPvalue <- min(tempRes$pvalue,na.rm=TRUE)
##   pvalues <- paste0(tempRes$pvalue,collapse=";")
##   nEnriched <- sum(tempRes$Enriched,na.rm=TRUE)
##   listofGene[[i]] <- data.frame(Gene=genes[i],meanLogFC,logFCs,minPvalue,pvalues,nEnriched)
## }
## geneTable <- do.call(rbind,listofGene)
## geneTable <- geneTable[order(geneTable$nEnriched,decreasing=TRUE),]


## ----include=FALSE------------------------------------------------------------
load("data/geneTable.RData")


## -----------------------------------------------------------------------------
geneTable[1:3,]


## ---- echo=TRUE,eval=FALSE,include=FALSE--------------------------------------
## myFQs <- c("/Users/thomascarroll/Downloads/PinAPL-py_demo_data/Control_R1_S14_L008_R1_001_x.fastq.gz", "/Users/thomascarroll/Downloads/PinAPL-py_demo_data/Control_R2_S15_L008_R1_001_x.fastq.gz", "/Users/thomascarroll/Downloads/PinAPL-py_demo_data/ToxA_R1_98_S2_L008_R1_001_x.fastq.gz", "/Users/thomascarroll/Downloads/PinAPL-py_demo_data/ToxA_R2_S2_L005_R1_001_x.fastq.gz", "/Users/thomascarroll/Downloads/PinAPL-py_demo_data/ToxB_R1_98_S4_L008_R1_001_x.fastq.gz", "/Users/thomascarroll/Downloads/PinAPL-py_demo_data/ToxB_R2_S4_L005_R1_001_x.fastq.gz")
## require(GenomicAlignments)
## counts <- list()
## counts2 <- list()
## stats <- list()
## 
## for(f in 1:length(myFQs)){
##   stats[[f]] <- align("GeCKO",myFQs[f],output_file = gsub(".fastq.gz",".bam",myFQs[f]),
##             nthreads=2,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1,maxMismatches = 0,indels = 0)
##   temp <- readGAlignments(gsub(".fastq.gz",".bam",myFQs[f]))
##   counts[[f]] <- data.frame(table(seqnames(temp[width(temp) == "20"])),row.names = "Var1")
##   counts2[[f]] <- data.frame(table(seqnames(temp)),row.names = "Var1")
## }
## myRes <- do.call(cbind,counts)
## 

