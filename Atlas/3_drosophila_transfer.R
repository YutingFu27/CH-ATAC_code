rm(list = ls())
library(data.table)
library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
#library(LSD)
library(reshape2)
addArchRThreads(threads = 12)
library(BSgenome.Dmelanogaster.UCSC.dm6)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Dmelanogaster.UCSC.dm6)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene, OrgDb = org.Dm.eg.db)

pbmc=readRDS('~/CH/CH-dm/res/pbmc_use_new.rds')


proj=readRDS('proj2_peak.rds')
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI_peak_res1",
  seRNA = pbmc,
  addToArrow = FALSE,
  groupRNA = "anno",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",force = TRUE
)

cM <- as.matrix(confusionMatrix(proj2$Clusters_peak_res1, proj2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
pal <- paletteDiscrete(values = pbmc$lineage)


pdf('transfer.pdf',height = 10,width = 10)
plotEmbedding(
  proj2, 
  colorBy = "cellColData", 
  embedding ="UMAP_peak_res1",
  name = "predictedGroup_Un", 
  pal = pal
)
dev.off()

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
saveRDS(proj2,file = './proj2_peak_intergrate.rds')