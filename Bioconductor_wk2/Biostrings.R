library(Biostrings)

dna1 = DNAString("ACGT-G")

dna2 = DNAStringSet(c("ACG","ACGT","ACGTT"))

IUPAC_CODE_MAP

dna1[2:4]
dna2[1:2]
dna2[[1]]

names(dna2) = paste0("seq",1:3)
dna2

sort(dna2)

rev(dna2)
rev(dna1)
reverseComplement(dna2)

translate(dna2)

alphabetFrequency(dna2)

letterFrequency(dna2, letter = "GC")

dinucleotideFrequency(dna2)

consensusMatrix(dna2)

#-----------------------------#

dnaseq <- DNAString("ACGTACGT")

vi = matchPattern(dnaseq,Scerevisiae$chrI)
ranges(vi)
shift(vi,10)

countPattern(dnaseq,Scerevisiae$chrI)

gr = vmatchPattern(dnaseq,Scerevisiae)
vi2 = Views(Scerevisiae, gr)

library(AnnotationHub)
ahub = AnnotationHub()
qh = query(ahub,c("sacCer2","genes"))
genes = qh[[1]]
prom = promoters(genes)
prom = trim(prom)
promViews = Views(Scerevisiaa,prom)
gcProm = letterFrequency(promViews,"GC",as.prob=TRUE)
plot(density(gcProm))

dnaseq == reverseComplement(dnaseq)

matchPWM()
pairwiseAlignment()
trimLRPatterns()


