library(AnnotationHub)
ah = AnnotationHub()
ah
ah[1]
ah[[1]]

ah = subset(ah, species = "Homo sapiens")
ah

query(ah, c("H3K4me3","Gm12878"))
ah2 = display(ah)
