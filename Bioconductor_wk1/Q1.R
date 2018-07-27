library(AnnotationHub)
library(GenomicRanges)
library(IRanges)
# Q1
ahub = AnnotationHub()
## Retrieve CpH Islands on genome "hg19"
gr = query(ahub,"CpG Islands")
CpG = gr[[1]]
autosome <- c(paste("chr", 1:22, sep=""))
length(keepSeqlevels(CpG,autosome,pruning.mode = "coarse"))
# Q2
length(keepSeqlevels(CpG,"chr4",pruning.mode = "coarse"))
# Q3
## Retrieve H3K4me3 histone modified H1 cell line on narrowPeak
qhs = query(ahub,c("H3K4me3","H1 cells"))
raw_H1_cell = qhs[[2]]
H1_cell = keepSeqlevels(raw_H1_cell,autosome,pruning.mode = "coarse")
sum(width(H1_cell))

# Q4
## Retrieve H3K27me3 histone modified H1 cell line on narrowPeak
qhs2 = query(ahub,c("H3K27me3","H1 cells"))
raw_H1_cell2 = qhs2[[2]]
H1_cell2 = keepSeqlevels(raw_H1_cell2,autosome,pruning.mode = "coarse")
mean(H1_cell2$signalValue)

# Q5
Bivalent = intersect(H1_cell,H1_cell2,ignore.strand = TRUE)
sum(width(Bivalent)) ## number of bases in both H3K4me3 and H3K27me3

# Q6
CpG_auto = keepSeqlevels(CpG,autosome,pruning.mode = "coarse")
length(subsetByOverlaps(Bivalent,CpG_auto,ignore.strand = TRUE))/length(Bivalent)

# Q7 
Bivalent_CpG = intersect(Bivalent,CpG_auto,ignore.strand = TRUE)
sum(width(Bivalent_CpG))/sum(width(CpG_auto))

# Q8 
CpG_resize = resize(CpG_auto, width(CpG_auto+10000), fix="center")
Bivalent_CpG_resize = intersect(Bivalent,CpG_resize,ignore.strand = TRUE)
sum(width(Bivalent_CpG_resize))

# Q9
sum(width(CpG_auto))/sum(as.numeric(seqlengths(CpG_auto)))

# Q10
## Bivalent on rows, CpG on cols
inOut = matrix(0, ncol=2,nrow=2)
colnames(inOut) = c("in","out")
rownames(inOut) = c("in","out")

inOut[1,1] = sum(width(intersect(Bivalent,CpG_auto,ignore.strand = TRUE)))

inOut[1,2] = sum(width(setdiff(Bivalent,CpG_auto,ignore.strand = TRUE))) ## number of bases only in Bivalent
inOut[2,1] = sum(width(setdiff(CpG_auto,Bivalent,ignore.strand = TRUE))) ## number of bases only in CpG

inOut[2,2] = sum(as.numeric(seqlengths(CpG_auto)))  - sum(inOut)
oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])