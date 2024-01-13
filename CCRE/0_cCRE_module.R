rm(list = ls())
setwd("~/CH/CH-zf/datu/bed/")


library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
#library(LSD)
library(reshape2)
addArchRThreads(threads =10)
library(BSgenome.Drerio.UCSC.danRer11)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Drerio.UCSC.danRer11)
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Drerio.UCSC.danRer11.refGene, OrgDb = org.Dr.eg.db)

load('~/CH/CH-zf/datu/Fig1/0_datu/proj2_datu_anno_pmat.rda')
proj2@peakSet

table(proj2$anno)

markersGS <- getMarkerFeatures(
  ArchRProj = proj2,
  useMatrix = "PeakMatrix",
  groupBy = 'anno',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.05")
save(markersGS,markerList,file = '~/CH/CH-zf/datu/Fig2/markerGS.rda')


heatmapPeaks <- ArchR::markerHeatmap(seMarker=markersGS,cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

pdf('~/CH/CH-zf/datu/Fig2/Goterm/dapeak_celltype_heatmap.pdf',height = 10,width = 15)
heatmapPeaks@ht_list$`Row Z-Scores
55729 features
PeakMatrix`
dev.off()


