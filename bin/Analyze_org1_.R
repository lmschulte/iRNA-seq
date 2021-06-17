#!/usr/bin/Rscript

# Variables from the bash script are stored in variable args
args <- commandArgs(trailingOnly=TRUE)

# Packages are loaded silently
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(splines))

# Parse arguments
Arguments <- as.data.frame(args)
Paired <- Arguments[1,1]
Dir <- Arguments[2,1]
Replication <- Arguments[3,1]
Disp <- Arguments[4,1]
Count <- Arguments[5,1]

AFiles <- data.frame()
BFiles <- data.frame()
Groups <- data.frame()
Counts <- data.frame()

i <- 6
while (Arguments[i,1] != "SEP") {
AFiles[i,1] <- Arguments[i,1]
i <- i + 1
}

i <- i + 1

while (Arguments[i,1] != "SEP2") {
BFiles[i,1] <- Arguments[i,1]
i <- i + 1
}

i <- i + 1

while (Arguments[i,1] != "SEP3") {
Groups[i,1] <- Arguments[i,1]
i <- i + 1
}

AFiles <- as.vector(AFiles[ !is.na(AFiles),])
BFiles <- as.vector(BFiles[ !is.na(BFiles),])
Groups <- as.vector(Groups[ !is.na(Groups),])
Paired <- as.vector(Paired)
Dir <- as.vector(Dir)
Replication <- as.vector(Replication)
Disp <- as.vector(Disp)
Count <- as.vector(Count)
Counts <- as.vector(Arguments[nrow(Arguments),1])

NoA <- length(AFiles)
NoB <- length(BFiles)
NoF <- NoA + NoB

# Import count data
cat("\t\t\t1.  Reading data and stitching genes together\n")
TMPDir <- paste(Dir, "/tmp/", sep="")
Counts <- paste(TMPDir, Counts, sep="")

if ( Count != "gene" && Count != "pol" ) {
Data <- read.delim(Counts)
Data$Substring <- substr( Data$Geneid, 5, nchar(as.character(Data$Geneid)))
Position <- data.frame(Position = regexpr("_", Data$Substring))
Data <- cbind(Data, Position)
Data$Gene <- substr( Data$Geneid, 1, Data$Position+3)
Data$Position <- regexpr("on", Data$Substring)
Data$Intron <- substr( Data$Substring, Data$Position+3, nchar(as.character(Data$Substring)))
Data <- Data[ !is.na(Data[,7]),]

# Calculate sums for each gene
Data_Sum <- aggregate( Data[,c(6:(7+NoF-1))], by=list(Data$Gene), FUN="sum")
colnames(Data_Sum)[1] <- "RefSeq"
rownames(Data_Sum) <- Data_Sum$RefSeq
Data_Sum$Sum <- rowSums(Data_Sum [,c(3:ncol(Data_Sum))])

} else if (Count == "gene" || Count == "pol") {
Data <- read.delim(Counts)

Data$Position <- regexpr("_", substr( Data[,1], 5, nchar(as.character(Data[,1]))))
Data[ Data$Position > -1, "Position"] <- Data[ Data$Position > -1, "Position"] + 3
Data[ Data$Position == "-1", "Position"] <- nchar(as.character(Data[Data$Position == "-1",1]))
Data[,1] <- substr(Data[,1], 1, Data$Position)

Data_Sum <- aggregate( Data[,c(6:(7+NoF-1))], by=list(Data[,1]), FUN="sum")
colnames(Data_Sum)[1] <- "RefSeq"
rownames(Data_Sum) <- Data_Sum$RefSeq
Data_Sum$Sum <- rowSums(Data_Sum [,c(3:ncol(Data_Sum))])
}

# Do differential expression analysis
cat("\t\t\t2.  Analyzing differential expression\n")
if (Paired == "yes") {
Cond <- factor(c(rep("A",NoA),rep("B",NoB)), levels=c("A","B"))
design = model.matrix(~Groups+Cond)
d <- DGEList(counts=data.matrix(Data_Sum [,c(3:(ncol(Data_Sum)-1))]), group=Cond)
d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)
glmfit <- glmFit(d, design)
Result <- glmLRT(glmfit)
Result <- topTags(Result, n=nrow(Data_Sum))
Result <- as.data.frame(Result)
colnames(Result) <- c("Log2_FC","logCPM","LR","Pval","Padj")
Result$RefSeq <- rownames(Result)
} else if (Paired == "no") {

if (Replication == "yes") {
Cond <- factor(c(rep("A",NoA),rep("B",NoB)), levels=c("A","B"))
d <- DGEList(counts=data.matrix(Data_Sum [,c(3:(ncol(Data_Sum)-1))]), group=Cond)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
Result <- exactTest(d)
Result <- topTags(Result, n=nrow(Data_Sum))
Result <- as.data.frame(Result)
colnames(Result) <- c("Log2_FC","logCPM","Pval","Padj")
Result$RefSeq <- rownames(Result)
} else if (Replication == "no") {
Cond <- factor(c(rep("A",NoA),rep("B",NoB)), levels=c("A","B"))
d <- DGEList(counts=data.matrix(Data_Sum [,c(3:(ncol(Data_Sum)-1))]), group=Cond)
d <- calcNormFactors(d)
d$common.dispersion <- as.numeric(as.character(Disp))
Result <- exactTest(d)
Result <- topTags(Result, n=nrow(Data_Sum))
Result <- as.data.frame(Result)
colnames(Result) <- c("Log2_FC","logCPM","Pval","Padj")
Result$RefSeq <- rownames(Result)
}
}

ANNFILE <- paste(Counts, ".ann", sep="")
Ann <- read.delim(ANNFILE, header=T)
Result <- merge(Result, Ann, by="RefSeq")
Result$txLength <- Result$End - Result$Start
Result <- merge(Result, Data_Sum[,c(1,2)], by="RefSeq")
colnames(Result)[ncol(Result)] <- "countLength"
CPM <- cpm(d,  normalized.lib.sizes=TRUE, log=T)
CPM <- as.data.frame(CPM)
CPM$RefSeq <- rownames(CPM)
Result <- merge(Result, CPM, by="RefSeq")
ResultOut <- paste(TMPDir, "out.txt", sep="")
#write.table(Result[,c(1,10,6,7,8,9,11,12,(13:(12+NoF)),2,3,4,5)], ResultOut, sep="\t", quote=F, col.names=T, row.names=F)
write.table(Result, ResultOut, sep="\t", quote=F, col.names=T, row.names=F)
