setwd("~/CH/CH-cross/cluster/monkey/tissue/")

library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
#library(LSD)
library(reshape2)
addArchRThreads(threads =1)
library(BSgenome.Mfascicularis.Ensemble.6.0)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mfascicularis.Ensemble.6.0)
genomeAnnotation$chromSizes

require(ChIPseeker)
library(GenomicFeatures)
load('../geneAnnotation.rda')

inputFiles <- list.files(pattern = 'bed.gz$')
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
proj<- ArchRProject(
  ArrowFiles = ArrowFiles,copyArrows = F,outputDirectory = '~/CH/CH-cross/cluster/monkey/pipeline/',geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation)

saveRDS(proj,file='~/CH/CH-cross/cluster/monkey/pipeline/proj.rds')

doubScores <- addDoubletScores(
  input = proj ,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
proj1 <- filterDoublets(doubScores) #23 of 150215
saveRDS(proj1,file='~/CH/CH-cross/cluster/monkey/pipeline/proj1.rds')

#anno
meta <- read.csv('./scATAC_meta.csv',sep = '\t')
head(meta)

meta$barcode=rownames(meta)
meta$barcode=colsplit(meta$barcode,'#',names = c('c1','c2'))$c2
meta[meta$Sample=='Liver_Normal',]$sample='Liver1'
meta[meta$Sample=='Liver',]$sample='Liver2'
table(meta$sample)

meta$bc=paste0(meta$sample,'#',meta$barcode)
head(meta$bc)

df <- as.data.frame(proj1@cellColData)
df$barcode=rownames(df)

bc.use=intersect(df$barcode,meta$bc)

proj2=proj1[proj1$cellNames%in%bc.use]

rownames(meta)=meta$bc
meta.use <- meta[proj2$cellNames,]
proj2$anno=meta.use$zhonghao
table(proj2$anno)

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

saveRDS(proj2,file = '~/CH/CH-cross/cluster/monkey/pipeline/proj2_bin.rds')

#peak-calling
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = 'Clusters_bin_res1',force = TRUE)
pathToMacs2 <- '/home/ggj/anaconda3/bin/macs2'
proj2<- addReproduciblePeakSet(
  ArchRProj = proj2, 
  groupBy = 'Clusters_bin_res1', 
  pathToMacs2=pathToMacs2,
  genomeSize = 2.7e9,excludeChr = 'MT'
)
proj2 <- addPeakMatrix(proj2)
saveRDS(proj2,file = '~/CH/CH-cross/cluster/monkey/pipeline/proj2_peakset.rds')
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

pdf('~/CH/CH-cross/cluster/monkey/pipeline/peak.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = addClusters_name, embedding = addUMAP_name)
dev.off()

saveRDS(proj2,file = '~/CH/CH-cross/cluster/monkey/pipeline/proj2_peak.rds')

proj2 <- addImputeWeights(proj2,reducedDims='IterativeLSI_peak_res1')


saveRDS(proj2,file = '~/CH/CH-cross/cluster/monkey/pipeline/proj2_impute.rds')

pdf('~/CH/CH-cross/cluster/monkey/pipeline/anno.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'anno', embedding = addUMAP_name)
dev.off()

pdf('~/CH/CH-cross/cluster/monkey/pipeline/sample.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = 'Sample', embedding = addUMAP_name)
dev.off()
