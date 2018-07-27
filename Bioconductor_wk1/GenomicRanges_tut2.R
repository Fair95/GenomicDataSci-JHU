library(GenomicRanges)
library(IRanges)

ir = IRanges(start = 1:3,width = 2)
df = DataFrame(ir = ir,score = rnorm(3))
df

df[1,1]
df$ir

df2 = data.frame(ir = ir)
df2

gr = GRanges(seqnames = c("chr1"),strand = c("+","-","+"),ranges = IRanges(start = c(1,3,5),width = 3))
values(gr) = DataFrame(score = rnorm(3))
gr
gr$score

df = data.frame(chr = "chr1",start = 1:3,end = 4:6,score =rnorm(3))
makeGRangesFromDataFrame(df,keep.extra.columns = TRUE)