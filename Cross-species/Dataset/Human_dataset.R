setwd("~/CH/CH-cross/cluster/Renbing/bed/")

library(Seurat)
library(ArchR)
library(reshape2)
addArchRThreads(threads = 12)
addArchRGenome("hg38")


sample.ls=gsub('.bed.gz','',inputFiles)


#bin.size 1000
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames=sample.ls,
  minTSS = 7,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  minFrags = 1000,
  TileMatParams = list(tileSize=1000)
)

ArrowFiles <- list.files(pattern = '*.arrow')
sample.ls=gsub('.arrow','',ArrowFiles )
proj<- ArchRProject(
  ArrowFiles = ArrowFiles,copyArrows = F,outputDirectory = '~/CH/CH-cross/cluster/Renbing/pipeline/')
saveRDS(proj,file='~/CH/CH-cross/cluster/Renbing/pipeline/proj.rds')

doubScores <- addDoubletScores(
  input = proj ,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
proj1 <- filterDoublets(doubScores) 
saveRDS(proj1,file = '~/CH/CH-cross/cluster/Renbing/pipeline/proj1.rds')

#add meta info
library(data.table)
meta <- fread('~/CH/CH-cross/cluster/Renbing/bed/all/meta/GSE184462_metadata.tsv.gz')
table(meta$`Life stage`)
meta=meta[meta$`Life stage`=='Adult',]
table(meta$cellID)


meta1 <- meta[meta$tissue%in%sample.ls,]
celltype <- table(meta1$`cell type`)%>%as.data.frame()

proj1$barcode=gsub('#','_1+',proj1$cellNames)
meta1=as.data.frame(meta1)
rownames(meta1)=meta1$cellID

bc_use=intersect(proj1$barcode,meta1$cellID)

proj2=proj1[proj1$barcode%in%bc_use]
proj2

meta1 <- meta1[proj2$barcode,]

meta1$cellID[1:5];proj2$barcode[1:5]

proj2$anno=meta1$`cell type`
table(proj2$anno)

sample20w=proj2$cellNames[sample(413540,200000)]
proj2=proj2[proj2$cellNames%in%sample20w]
table(proj2$Sample)%>%as.data.frame()->df
summary(df$Freq)
table(proj2$anno)
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

saveRDS(proj2,file = '~/CH/CH-cross/cluster/Renbing/pipeline/proj2_bin.rds')

#peak-calling
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = 'Clusters_bin_res1',force = TRUE)
pathToMacs2 <- '/home/ggj/anaconda3/bin/macs2'
proj2<- addReproduciblePeakSet(
  ArchRProj = proj2, 
  groupBy = 'Clusters_bin_res1', 
  pathToMacs2=pathToMacs2,
  genomeSize = 2.7e9
)
saveRDS(proj2,file = '~/CH/CH-cross/cluster/Renbing/pipeline/proj2_peakset.rds')
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

pdf('~/CH/CH-cross/cluster/Renbing/pipeline/peak.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = addClusters_name, embedding = addUMAP_name)
dev.off()

pdf('~/CH/CH-cross/cluster/Renbing/pipeline/anno.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'anno', embedding = addUMAP_name)
dev.off()

saveRDS(proj2,file = '~/CH/CH-cross/cluster/Renbing/pipeline/proj2_peak.rds')

proj2 <- addImputeWeights(proj2,reducedDims='IterativeLSI_peak_res1')


saveRDS(proj2,file = '~/CH/CH-cross/cluster/Renbing/pipeline/proj2_impute.rds')

proj2$sample=reshape2::colsplit(proj2$Sample,'_',names = c('c1','c2'))$c1
table(proj2$sample)

pdf('~/CH/CH-cross/cluster/Renbing/pipeline/tissue.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'sample', embedding = addUMAP_name)
dev.off()
