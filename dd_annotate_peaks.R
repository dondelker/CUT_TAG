# Load libraries
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPpeakAnno)
library(biomaRt)

# Add genome information
genome.hg38.ensembl = BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(genome.hg38.ensembl) = "Ensembl"
data("TSS.human.GRCh38")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Read in Granges list of CUT&Tag peak regions
CTpeaks = GRanges(read.csv("FB_precision_coordinates.csv"))

# Annotate peak regions
annoPeaks = annotatePeakInBatch(CTpeaks, AnnotationData = TSS.human.GRCh38)
annoPeaks = addGeneIDs(annoPeaks, mart=mart, IDs2Add=c("hgnc_symbol"))

# Retrieve DNA sequence for each peak region
annotate.gene.loc = function(annoPeaks) {
  seqlevelsStyle(annoPeaks) = "Ensembl"
  annoPeaks$seq = as.character(getSeq(genome.hg38.ensembl, annoPeaks))
}
annoPeaks = sort(annoPeaks)
add_sequence = annotate.gene.loc(annoPeaks)

# Combine dataframes
annoPeaks <- as.data.frame(annoPeaks)
sequence <- as.data.frame(add_sequence)
geneFinal <- cbind(annoPeaks, sequence)

# Write to file
write.table(geneFinal, "FB_precision_1kb.txt", sep = "\t")