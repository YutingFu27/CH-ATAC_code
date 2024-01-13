setwd('~/CH/CH-ew/bed_total/')


library(Seurat)
library(ArchR)
load('~/CH/CH-ew/res/ew_archR.rda')

library(data.table)
addArchRThreads(threads =18)


pbmc <- readRDS('~/CH/CH-ew/RNA/ew_use.rds')
pbmc
table(pbmc$cell_type)

proj2 <- addGeneIntegrationMatrix(
  ArchRProj = proj2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI_peak_res1",
  seRNA = pbmc,
  addToArrow = FALSE,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",force = TRUE
)

cM <- as.matrix(confusionMatrix(proj2$Clusters_peak_res1, proj2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
pal <- paletteDiscrete(values = pbmc$cell_type)

pdf('transfer.pdf',height = 10,width = 10)
plotEmbedding(
  proj2, 
  colorBy = "cellColData", 
  embedding ="UMAP_peak_res1",
  name = "predictedGroup_Un", 
  pal = pal
)
dev.off()

pdf('pred_score.pdf',height = 10,width = 10)
hist(proj2$predictedScore_Un,xlab = 'prediction score',col = 'light blue',xlim = c(0,1),main = 'CHATAC-zebrafish')
abline(v=0.5,col='red',lwd=2,lty=2)
dev.off()

#pheatmap
cM <- cM/Matrix::rowSums(cM)
pheatmap::pheatmap(
  mat=as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = 'black',cluster_rows = F,cluster_cols = F
)