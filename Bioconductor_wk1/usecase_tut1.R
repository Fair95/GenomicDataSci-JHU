library(AnnotationHub)

ahub = AnnotationHub()
ahub = subset(ahub, species == "Homo sapiens")
qhs = query(ahub, c("H3K4me3","Gm12878"))
gr1 = qhs[[2]]
gr2 = qhs[[4]]
summary(width(gr1))
summary(width(gr2))
table(width(gr2))

## Retrive peaks
peaks = gr2
qhs[4]

# Find reference genes
qhs = query(ahub,"RefSeq")
qhs$genomes
genes = qhs[[1]]

# Find promoters in ref genes
prom = promoters(genes)
table(width(prom))

promoters

peaks

## check whether the peaks enrich in promoters
ov = findOverlaps(prom,peaks)
length(unique(queryHits(ov))) ## number of unique hits

subsetByOverlaps(peaks,prom,ignore.strand = TRUE)/length(peaks) ## percentage of mapping peaks

subsetByOverlaps(prom,peaks,ignore.strand = TRUE)/length(prom) ## percentage of mapping prom

sum(width(reduce(peaks, ignore.strand = TRUE)))/10^6

sum(width(reduce(prom, ignore.strand = TRUE)))/10^6

sum(width(intersect(peaks,prom,ignore.strand = TRUE)))/10^6

## peaks on rows, prom on cols
inOut = matrix(0, ncol=2,nrow=2)
colnames(inOut) = c("in","out")
rownames(inOut) = c("in","out")

inOut[1,1] = sum(width(intersect(peaks,prom,ignore.strand = TRUE))) ## number of bases in both peaks and prom
inOut[1,2] = sum(width(setdiff(peaks,prom,ignore.strand = TRUE))) ## number of bases only in peaks
inOut[2,1] = sum(width(setdiff(prom,peaks,ignore.strand = TRUE))) ## number of bases only in prom

colSums(inOut) ## =sum(width(reduce(prom, ignore.strand = TRUE)))
rowSums(inOut) ## =sum(width(reduce(peaks, ignore.strand = TRUE)))

inOut[2,2] = 3*10^9 - sum(inOut) ## bases not in either of peaks or prom
inOut
#fisher.test(inOut)$statistic will report error because human genome 3*10^9 >.Machine$integer.max
oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
oddsRatio 

## not all human genes are mappable
## lets assume only half mapped
inOut[2,2] = 1.5 *10^9 - sum(inOut)

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
oddsRatio