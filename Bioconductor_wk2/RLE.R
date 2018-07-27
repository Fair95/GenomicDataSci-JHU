library(GenomicRanges)

rl = Rle(c(1,1,1,1,1,1,2,2,2,2,2,4,4,2))

runLength(rl)
runValue(rl)

as.numeric(rl)

ir = IRanges(start = c(2,8),width = 4)
aggregate(rl,ir,FUN = mean)
vec = as.numeric(rl)
mean(vec[2:5])
mean(vec[8:11])

ir = IRanges(start=1:5,width = 3)
coverage(ir)

rl
slice(rl,2)
slice(rl,3)

gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:4,width = 3))
gr2 <- GRanges(seqnames = "chr2", ranges = IRanges(start = 1:4,width = 3))
gL = GRangesList(gr1 = gr1, gr2 = gr2)
gL
gL[1]
gL[[1]]
gL$gr1
start(gL)
seqnames(gL)

elementLengths(gL)
sapply(gL,length)
shift(gl, 10)

findOverlaps(gL, gr2)