# Variables from the bash script are stored in variable args
args <- commandArgs(trailingOnly=TRUE)

if (args[2] == "1") {
DumpFile <- paste(args[1],"/Gene.Dump",sep="")
Dump <- read.delim(DumpFile, dec=",", header=F)

# Collapse gene lists based on Annotation
Dump$Test <- regexpr("_", Dump$V3)
Dump <- Dump[ Dump$Test == "-1",c(1:(ncol(Dump)-1))]
Dump2 <- Dump
Dump$Length <- Dump$V6-Dump$V5
Dump$Sub <- substr(Dump$V2,0,2)
Dump <- Dump[ Dump$Sub %in% c("AT"),]
Dump2 <- Dump2[ !(Dump2$V13 %in% Dump$V13),]
Dump2$Length <- Dump2$V6-Dump2$V5
Dump2$Sub <- substr(Dump2$V2,0,2)
Dump <- rbind(Dump, Dump2)

Dump <- Dump[ order(Dump$V13, -Dump$Length),]
Dump <- Dump[ duplicated(Dump$V13)==F,]
Dump <- Dump[ order(Dump$V2, -Dump$Length),]
Dump <- Dump[ duplicated(Dump$V2)==F,]

# Make filename for output
DumpFile <- paste(args[1],"/Gene.Dump.Inc",sep="")

# Output regions in included genes
write.table( Dump, DumpFile, sep="\t", quote=F, col.names=F, row.names=F)
}

if (args[2] == "2") {
DumpFile <- paste(args[1],"/Gene.Dump.Inc",sep="")
Dump <- read.delim(DumpFile, dec=",", header=F)

RegionFile <- paste(args[1],"/Intron.filt.unstranded.bed",sep="")
Region <- read.delim(RegionFile, dec=",", header=F)
Region <- merge(Region, Dump[,c(2,13,18)], all.x=T, by.x="V4", by.y="V13")
Region$Count <- 1
Count <- aggregate(Region$Count, by=list(Region$V2.y), FUN="sum")
Region <- merge(Region, Count, by.x="V2.y", by.y="Group.1", all.x=T)
write.table(Region[ Region$V18 %in% c("AT"), c(1,3,4,5,7,10)], RegionFile, quote=F, row.names=F, col.names=F, sep="\t")

RegionFile <- paste(args[1],"/Intron.filt.stranded.bed",sep="")
Region <- read.delim(RegionFile, dec=",", header=F)
Region <- merge(Region, Dump[,c(2,13,18)], all.x=T, by.x="V4", by.y="V13")
Region$Count <- 1
Count <- aggregate(Region$Count, by=list(Region$V2.y), FUN="sum")
Region <- merge(Region, Count, by.x="V2.y", by.y="Group.1", all.x=T)
write.table(Region[ Region$V18 %in% c("AT"), c(1,3,4,5,7,10)], RegionFile, quote=F, row.names=F, col.names=F, sep="\t")

RegionFile <- paste(args[1],"/Exon.filt.unstranded.bed",sep="")
Region <- read.delim(RegionFile, dec=",", header=F)
Region <- merge(Region, Dump[,c(2,13,18)], all.x=T, by.x="V4", by.y="V13")
Region$Count <- 1
Count <- aggregate(Region$Count, by=list(Region$V2.y), FUN="sum")
Region <- merge(Region, Count, by.x="V2.y", by.y="Group.1", all.x=T)
write.table(Region[ Region$V18 %in% c("AT"), c(1,3,4,5,7,10)], RegionFile, quote=F, row.names=F, col.names=F, sep="\t")

RegionFile <- paste(args[1],"/Exon.filt.stranded.bed",sep="")
Region <- read.delim(RegionFile, dec=",", header=F)
Region <- merge(Region, Dump[,c(2,13,18)], all.x=T, by.x="V4", by.y="V13")
Region$Count <- 1
Count <- aggregate(Region$Count, by=list(Region$V2.y), FUN="sum")
Region <- merge(Region, Count, by.x="V2.y", by.y="Group.1", all.x=T)
write.table(Region[ Region$V18 %in% c("AT"), c(1,3,4,5,7,10)], RegionFile, quote=F, row.names=F, col.names=F, sep="\t")

RegionFile <- paste(args[1],"/Gene.filt.unstranded.bed",sep="")
Region <- read.delim(RegionFile, dec=",", header=F)
Region <- merge(Region, Dump[,c(2,13,18)], all.x=T, by.x="V4", by.y="V13")
Region$Count <- 1
Count <- aggregate(Region$Count, by=list(Region$V2.y), FUN="sum")
Region <- merge(Region, Count, by.x="V2.y", by.y="Group.1", all.x=T)
write.table(Region[ Region$V18 %in% c("AT"), c(1,3,4,5,7,10)], RegionFile, quote=F, row.names=F, col.names=F, sep="\t")

RegionFile <- paste(args[1],"/Gene.filt.stranded.bed",sep="")
Region <- read.delim(RegionFile, dec=",", header=F)
Region <- merge(Region, Dump[,c(2,13,18)], all.x=T, by.x="V4", by.y="V13")
Region$Count <- 1
Count <- aggregate(Region$Count, by=list(Region$V2.y), FUN="sum")
Region <- merge(Region, Count, by.x="V2.y", by.y="Group.1", all.x=T)
write.table(Region[ Region$V18 %in% c("AT"), c(1,3,4,5,7,10)], RegionFile, quote=F, row.names=F, col.names=F, sep="\t")
}
