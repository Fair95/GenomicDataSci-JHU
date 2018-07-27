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
qhs = query(ahub,c("H3K4me3","H1 cells"))
raw_H1_cell = qhs[[2]]
H1_cell_chr22 = keepSeqlevels(raw_H1_cell,"chr22",pruning.mode = "coarse")
h1view = View(Hsapiens,ranges(H1_cell_chr22))
param = new("BSParams", X = h1_view, FUN = letterFrequency)
mean(unlist(bsapply(param,"GC",,as.prob = TRUE)))

# Q3
cor()
