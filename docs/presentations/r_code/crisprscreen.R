params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Crispr Screenining analysis in R

---
"    
  )
  
}



## ----downloadFile-------------------------------------------------------------
download.file("https://github.com/LewisLabUCSD/PinAPL-Py/archive/master.zip","pieappleData.zip")
unzip("pieappleData.zip")


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



## ----fa1----------------------------------------------------------------------
download.file("http://pinapl-py.ucsd.edu/example-data","TestData.zip")
unzip("TestData.zip")


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

temp <- readGAlignments(gsub(".fastq.gz",".bam",myFQs))
temp


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
table(grepl("29S20M1S",temp@cigar))
temp <- temp[grepl("29S20M1S",temp@cigar)]

# temp data.frame(table(seqnames(temp[width(temp) == "20"])),row.names = "Var1")


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
counts <- data.frame(table(seqnames(temp)),row.names = "Var1")
counts <- data.frame(table(seqnames(temp[width(temp) == "20"])),row.names = "Var1")
counts[1:2,,drop=FALSE]


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
download.file("http://pinapl-py.ucsd.edu/run/download/example-run","Results.zip")
unzip("Results.zip")


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
ss <- read.delim("New_Example_1580342908_4/Analysis/01_Alignment_Results/ReadCounts_per_sgRNA/Control_1_GuideCounts.txt")
new <- merge(ss,counts,by.x=1,by.y=0)
plot(new$X0,new$Freq)

