library(GenomicRanges)
gr = GRanges(seqnames = c("chr1","chr2","chr2"),strand = c("+","-","+"),ranges = IRanges(start = c(1,3,5),width = 3))

## keep only chr1
seqlevels(gr,pruning.mode = "coarse") = "chr1"

gr = GRanges(seqnames = c("chr1","chr2","chr2"),strand = c("+","-","+"),ranges = IRanges(start = c(1,3,5),width = 3))
dropSeqlevels(gr,"chr2",pruning.mode = "coarse")

keepSeqlevels(gr,"chr1",pruning.mode = "coarse")

## convert and unify styles
newStyle = mapSeqlevels(seqlevels(gr),"NCBI")
newStyle

gr = renameSeqlevels(gr, newStyle)