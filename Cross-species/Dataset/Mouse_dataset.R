rm(list = ls())
setwd("~/sci-mouseATAC/sci1/")

library(Seurat)
library(ArchR)
library(reshape2)
addArchRThreads(threads = 8)
addArchRGenome("mm9")

ArrowFiles <- list.files(pattern = '*.arrow')
ArrowFiles=ArrowFiles[-1]
proj<- ArchRProject(
  ArrowFiles = ArrowFiles,copyArrows = F,outputDirectory = '~/sci-mouseATAC/sci1/pipeline/')

saveRDS(proj,file='~/sci-mouseATAC/sci1/pipeline/proj.rds')

doubScores <- addDoubletScores(
  input = proj ,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
proj1 <- filterDoublets(doubScores) #23 of 150215
saveRDS(proj1,file='~/sci-mouseATAC/sci1/pipeline/proj1.rds')
proj1<-readRDS('~/sci-mouseATAC/sci1/pipeline/proj1.rds')

meta <- fread('./cell_metadata.tissue_freq_filtered.txt')
meta=as.data.frame(meta)
meta$bc=paste0(meta$tissue.replicate,'#',meta$cell)
rownames(meta)=meta$bc

use.bc=intersect(rownames(meta),proj1$cellNames)
proj2=proj1[proj1$cellNames%in%use.bc]

meta.use=meta[proj2$cellNames,]
table(meta.use$cell_label)
meta.use$anno=meta.use$cell_label

proj2$anno=meta.use$anno
table(proj2$anno)
proj2=proj2[proj2$anno%ni%c('Unknown','Collisions')]

proj2
##reduce dimension
res=1
addIterativeLSI_name=paste0('IterativeLSI_bin_res',res)
proj2 <- addIterativeLSI(
  ArchRProj = proj2,
  useMatrix = "TileMatrix", 
  name = addIterativeLSI_name, 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = res, 
    sampleCells = NULL, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:50,force = TRUE
)
# ##clustering
addClusters_name=paste0('Clusters_bin_res',res)
proj2 <- addClusters(
  input = proj2,
  reducedDims = addIterativeLSI_name,
  method = "Seurat",
  name = addClusters_name,
  resolution = res,force = TRUE
)
##umap
addUMAP_name=paste0('UMAP_bin_res',res)
proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = addIterativeLSI_name, 
  name =addUMAP_name, 
  nNeighbors = 30, 
  minDist = 0.8, 
  metric = "cosine",force = TRUE
)
saveRDS(proj2,file = '~/sci-mouseATAC/sci1/pipeline/proj2_bin.rds')

proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = 'Clusters_bin_res1',force = TRUE)
pathToMacs2 <- '/home/ggj/anaconda3/bin/macs2'
proj2<- addReproduciblePeakSet(
  ArchRProj = proj2, 
  groupBy = 'Clusters_bin_res1', 
  pathToMacs2=pathToMacs2,
  genomeSize = 1.87e9
)
proj2 <- addPeakMatrix(proj2)

#peak-cell cluster
res=1
addIterativeLSI_name=paste0('IterativeLSI_peak_res',res)
proj2 <- addIterativeLSI(
  ArchRProj = proj2,
  useMatrix = "PeakMatrix", 
  name = addIterativeLSI_name, 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = res, 
    sampleCells = NULL, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:50
)
##clustering
addClusters_name=paste0('Clusters_peak_res',res)
proj2 <- addClusters(
  input = proj2,
  reducedDims = addIterativeLSI_name,
  method = "Seurat",
  name = addClusters_name,
  resolution = res,
  force = TRUE,maxClusters = 100
)
##umap
addUMAP_name=paste0('UMAP_peak_res',res)
proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = addIterativeLSI_name, 
  name = addUMAP_name, 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",force = TRUE
)

pdf('~/sci-mouseATAC/sci1/pipeline/peak.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = addClusters_name, embedding = addUMAP_name)
dev.off()

pdf('~/sci-mouseATAC/sci1/pipeline/anno.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'anno', embedding = addUMAP_name)
dev.off()
saveRDS(proj2,file = '~/sci-mouseATAC/sci1/pipeline/proj2_peak.rds')

proj2 <- addImputeWeights(proj2,reducedDims='IterativeLSI_peak_res1')
proj2$sample=reshape2::colsplit(proj2$Sample,'_',names = c('c1','c2'))$c1
table(proj2$sample)
saveRDS(proj2,file = '~/sci-mouseATAC/sci1/pipeline/proj2_impute.rds')

pdf('~/sci-mouseATAC/sci1/pipeline/tissue.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'sample', embedding = addUMAP_name)
dev.off()

setwd('~/sci-mouseATAC/sci1/')
proj2 <- readRDS('~/sci-mouseATAC/sci1/pipeline/1_subcluster/proj2_anno_all.rds')
pdf('~/sci-mouseATAC/sci1/pipeline/mm9_anno_our.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'cell_type', embedding = addUMAP_name)
dev.off()
