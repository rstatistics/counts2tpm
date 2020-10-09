library(GenomicFeatures)

gtffile <- "gencode.v22.annotation.gtf"

txdb <- makeTxDbFromGFF(gtffile, format="gtf")

ebg <- exonsBy(txdb, by="gene")

Counts <- read.table(file = "Counts.txt", header = TRUE, row.names = 1)

ebgList <- sum(width(reduce(ebg)))

genes <- intersect(rownames(Counts), names(ebgList))

Length <- as.vector(ebgList[genes])

TPM <- t(t(Counts / t(Length)) * 1e6 / colSums(Counts / t(Length)))

FPKM <- t(t(Counts / t(Length)) * 1e9 / colSums(Counts))

write.table(TPM, file="TPM.txt", append=FALSE, quote=FALSE, sep="\t", row.names=1, col.names=1)

write.table(FPKM, file="FPKM.txt", append=FALSE, quote=FALSE, sep="\t", row.names=1, col.names=1)
