rm(list = ls())
setwd("~/CH/CH-zf/datu/bed/")

library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
#library(LSD)
library(reshape2)
addArchRThreads(threads =23)
library(BSgenome.Drerio.UCSC.danRer11)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Drerio.UCSC.danRer11)
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Drerio.UCSC.danRer11.refGene, OrgDb = org.Dr.eg.db)

proj=readRDS('proj2_peak.rds')
load('~/CH/CH-zf/datu/bed/anno/zclrnasample.rdata')

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI_bin_res1",
  seRNA = pbmc,
  addToArrow = FALSE,
  groupRNA = "anno",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(proj$Clusters_bin_res1, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
pal <- paletteDiscrete(values = pbmc$anno)

pdf('transfer.pdf',height = 10,width = 10)
plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  embedding ="UMAP_bin_res1",
  name = "predictedGroup_Un", 
  pal = pal
)
dev.off()

proj2 <- addGeneIntegrationMatrix(
  ArchRProj = proj2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI_peak_res1",
  seRNA = pbmc,
  addToArrow = FALSE,
  groupRNA = "anno",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
saveRDS(proj2,file = './proj2_peak_intergrate.rds')
cM <- as.matrix(confusionMatrix(proj$Clusters_peak_res1, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
pal <- paletteDiscrete(values = pbmc$anno)

#hist
pdf('pred_score.pdf',height = 10,width = 10)
hist(proj2$predictedScore_Un,xlab = 'prediction score',col = 'light blue',xlim = c(0,1),main = 'CHATAC-earthworm')
abline(v=0.5,col='red',lwd=2,lty=2)
dev.off()

#pheatmap
cM <- cM/Matrix::rowSums(cM)
pheatmap::pheatmap(
  mat=as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = 'black',cluster_rows = F,cluster_cols = F
)