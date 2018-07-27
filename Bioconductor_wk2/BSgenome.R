library(BSgenome)
available.genomes()


Scerevisiae = available.genomes()[[77]]
Scerevisiae

param = new("BSParams", X = Scerevisiae, FUN = letterFrequency)
sum(unlist(bsapply(param,"GC")))/sum(seqlengths(Scerevisiae))

unlist(bsapply(param,"GC",as.prob = TRUE))
