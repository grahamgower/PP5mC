args = commandArgs(trailingOnly=TRUE)

fileN<-args[1]

#fileN<-"A16121-MCS34-2.Cow_UMD3_1.pileOmeth.CpG.binStats"

dat<-as.data.frame(read.table(fileN, header=F, na.string="."))

dhis<-hist(dat[,13], breaks=60, xlim=c(0,30))
chrom<-unique(dat[,1])

outFile<-gsub("[.]binStats",".pdf",fileN)

pdf(outFile)
hist(dat[which(dat[,13]>=5),4], xlab="Methylation level", col="gray", ylab="Counts", main="Global methylation level", sub="Cs covered by >= 5 reads")
hist(dat[which(dat[,13]>=3),4], xlab="Methylation level", ylab="Counts", col="gray", main="Global methylation level", sub="Cs covered by >= 3 reads")
hist(dat[,4], xlab="Methylation level", ylab="Counts", col="gray", main="Global methylation level", sub="All Cs")

plot(dhis$counts, log="y", type="h", xlab="Coverage level per 1kb window", ylab="Counts", main="Global coverage")

par(mfrow=c(3,3))
for(i in 1:length(chrom)){
	plot(dat[which(dat[,1]==chrom[i]),2], dat[which(dat[,1]==chrom[i]),13], type="l",xlab=paste(chrom[i],", 1kb bins"), ylab="Coverage level (log scale)", col="deepskyblue4", log="y")
}

par(mfrow=c(2,1))
boxplot(V4~V1, data=dat, notch=TRUE, las=2, ylab="Methylation level", main="Methylation level per chromosome", col="gold")
boxplot(V13~V1, data=dat, log="y", las=2, ylab="Coverage level (log scale)", main="Coverage level per chromosome", col="gold")

dev.off()
