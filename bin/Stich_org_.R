#!/usr/bin/Rscript

# Variables from the bash script are stored in variable args
args <- commandArgs(trailingOnly=TRUE)

# Parse arguments
Arguments <- as.data.frame(args)
Dir <- Arguments[1,1]
Count <- Arguments[2,1]
String <- Arguments[3,1]
Norm <- Arguments[4,1]

AFiles <- data.frame()

i <- 5
while (i <= nrow(Arguments)) {
AFiles[i,1] <- Arguments[i,1]
i <- i + 1
}

AFiles <- as.vector(AFiles[ !is.na(AFiles),])
Dir <- as.vector(Dir)
Count <- as.vector(Count)
String <- as.vector(String)
NoA <- length(AFiles)
Norm <- as.vector(Norm)

# Import count data
TMPDir <- paste(Dir, "/tmp/", sep="")
Counts <- paste(TMPDir, String, sep="")

if ( Count != "gene" && Count != "pol" ) {
Data <- read.delim(Counts)
Data$Substring <- substr( Data$Geneid, 5, nchar(as.character(Data$Geneid)))
Position <- data.frame(Position = regexpr("_", Data$Substring))
Data <- cbind(Data, Position)
Data$Gene <- substr( Data$Geneid, 1, Data$Position+3)
Data$Position <- regexpr("on", Data$Substring)
Data$Intron <- substr( Data$Substring, Data$Position+3, nchar(as.character(Data$Substring)))
Data <- Data[ !is.na(Data[,7]),]

# Calculate sums for each gene and remove genes without tags
Data_Sum <- aggregate( Data[,c(6:(7+NoA-1))], by=list(Data$Gene), FUN="sum")
colnames(Data_Sum)[1] <- "RefSeq"
rownames(Data_Sum) <- Data_Sum$RefSeq
Data_Sum$Sum <- rowSums(Data_Sum [,c(3:ncol(Data_Sum))])

} else if (Count == "gene" || Count == "pol") {
Data <- read.delim(Counts)

Data$Position <- regexpr("_", substr( Data[,1], 5, nchar(as.character(Data[,1]))))
Data[ Data$Position > -1, "Position"] <- Data[ Data$Position > -1, "Position"] + 3
Data[ Data$Position == "-1", "Position"] <- nchar(as.character(Data[Data$Position == "-1",1]))
Data[,1] <- substr(Data[,1], 1, Data$Position)

Data_Sum <- aggregate( Data[,c(6:(7+NoA-1))], by=list(Data[,1]), FUN="sum")
colnames(Data_Sum)[1] <- "RefSeq"
rownames(Data_Sum) <- Data_Sum$RefSeq
Data_Sum$Sum <- rowSums(Data_Sum [,c(3:ncol(Data_Sum))])
}

if ( Norm == "no" ) {
ANNFILE <- paste(Counts, ".ann", sep="")
Ann <- read.delim(ANNFILE, header=T)
Data_Sum <- merge(Data_Sum, Ann, by="RefSeq")
Data_Sum$txLength <- Data_Sum$End - Data_Sum$Start
colnames(Data_Sum)[2] <- "countLength"
ResultOut <- paste(TMPDir, "out.txt", sep="")
write.table(Data_Sum[,c(1,(2+NoA+6),(2+NoA+2),(2+NoA+4),(2+NoA+5),(2+NoA+3),(2+NoA+7),2,3:(2+NoA))], ResultOut, sep="\t", quote=F, col.names=T, row.names=F)
}

if ( Norm == "yes" ) {
suppressPackageStartupMessages(require(edgeR))

d <- DGEList(counts=data.matrix(Data_Sum[,c(3:(2+NoA))]))
d <- calcNormFactors(d)
CPM <- cpm(d,  normalized.lib.sizes=TRUE, log=T)

CPM <- as.data.frame(CPM)
CPM$RefSeq <- rownames(CPM)
Data_Sum <- merge(Data_Sum[,c(1,2)], CPM, by="RefSeq")

ANNFILE <- paste(Counts, ".ann", sep="")
Ann <- read.delim(ANNFILE, header=T)
Data_Sum <- merge(Data_Sum, Ann, by="RefSeq")
Data_Sum$txLength <- Data_Sum$End - Data_Sum$Start
colnames(Data_Sum)[2] <- "countLength"
ResultOut <- paste(TMPDir, "out.txt", sep="")
write.table(Data_Sum[,c(1,(2+NoA+5),(2+NoA+1),(2+NoA+3),(2+NoA+4),(2+NoA+2),(2+NoA+6),2,3:(2+NoA))], ResultOut, sep="\t", quote=F, col.names=T, row.names=F)
}
