library(GenomicRanges)

gr = GRanges(seqnames = c("chr1"),strand = c("+","-","+"),ranges = IRanges(start = c(1,3,5),width = 3))

flank(gr,5)

promoters(gr)

seqinfo(gr)

seqlengths(gr) = c("chr1" = 10)

seqinfo(gr)

gaps(gr)

seqlevels(gr) = c("chr1","chr2")

seqnames(gr) = c("chr1","chr2","chr1")


sort(gr)
seqlevels(gr) = c("chr2","chr1")
sort(gr)

genome(gr) = "hg19"
seqinfo(gr)
