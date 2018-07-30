library(AnnotationHub)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
# Q1 
GC <- letterFrequency(Hsapiens$chr22, "GC")
ACGT <- letterFrequency(Hsapiens$chr22,"ACGT")
ratio = GC/ACGT

# Q2
ahub = AnnotationHub()
qhs = query(ahub,c("H3K27me3","H1 cells"))
raw_H1_cell = qhs[[2]]
H1_cell_chr22 = keepSeqlevels(raw_H1_cell,"chr22",pruning.mode = "coarse")
h1_view = Views(Hsapiens,H1_cell_chr22)
#param = new("BSParams", X = h1_view, FUN = letterFrequency)
#mean(unlist(bsapply(param,"GC",,as.prob = TRUE)))
GC <- letterFrequency(h1_view,"GC",as.prob = TRUE)
mean(GC)
# Q3
signal = matrix(H1_cell_chr22$signalValue)
cor(signal,GC)

# Q4
qhs = query(ahub,c("H3K27me3","H1 cells","fc.signal"))
fc = qhs[[1]]
fc_chr22 = keepSeqlevels(fc,"chr22",pruning.mode = "coarse")
