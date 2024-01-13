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

#loading fragment files
file.list=list.files(pattern = '_zf_sort.bed.gz$')
sample.list=gsub('_zf_sort.bed.gz','',file.list)

ArrowFiles <- createArrowFiles(
  inputFiles =file.list,
  sampleNames=sample.list,
  minTSS = 7,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  minFrags = 1000,
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation
)

ArrowFiles<-list.files(pattern = '*.arrow')
proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation,copyArrows = F)


##rm doublets
doubScores <- addDoubletScores(
  input = proj1 ,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
proj1 <- filterDoublets(doubScores) 

#QC plot
df <- getCellColData(proj1, select = c("log10(nFrags)", "TSSEnrichment"))
ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment")


##cell-bin clustering
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


addClusters_name=paste0('Clusters_bin_res',res)
proj2 <- addClusters(
  input = proj2,
  reducedDims = addIterativeLSI_name,
  method = "Seurat",
  name = addClusters_name,
  resolution = res,force = TRUE
)

addUMAP_name=paste0('UMAP_bin_res',res)
proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = addIterativeLSI_name, 
  name =addUMAP_name, 
  nNeighbors = 30, 
  minDist = 0.8, 
  metric = "cosine",force = TRUE
)

pdf('bin.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = addClusters_name, embedding = addUMAP_name)
dev.off()

saveRDS(proj2,file = './proj2_bin.rds')

#peak-calling
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = 'Clusters_bin_res1',force = TRUE)
pathToMacs2 <- '/home/ggj/anaconda3/envs/SAMap/bin/macs2'
proj2<- addReproduciblePeakSet(
  ArchRProj = proj2, 
  groupBy = 'Clusters_bin_res1', 
  pathToMacs2=pathToMacs2,
  genomeSize=1.5e9 #zf genome size
)
proj2 <- addPeakMatrix(proj2)

saveRDS(proj2,file = './proj2_pmat.rds')


#cell-peak clustering
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

addClusters_name=paste0('Clusters_peak_res',res)
proj2 <- addClusters(
  input = proj2,
  reducedDims = addIterativeLSI_name,
  method = "Seurat",
  name = addClusters_name,
  resolution = res,
  force = TRUE,maxClusters = 100
)

addUMAP_name=paste0('UMAP_peak_res',res)
proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = addIterativeLSI_name, 
  name = addUMAP_name, 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",force = TRUE
)
saveRDS(proj2,file = 'proj2_peak.rds')
pdf('peak.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = addClusters_name, embedding = addUMAP_name)
dev.off()

#marker gene list
markersGS <- getMarkerFeatures(
  ArchRProj = proj2,
  useMatrix = "GeneScoreMatrix",
  groupBy = addClusters_name,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")

## output markerlist for clusters
marker=markerList@listData
names(marker)=c(1:length(markerList))
use<-data.frame(matrix(ncol = 3,nrow = 0))
colnames(use)=c('cluster','gene','Log2FC')
for (i in 1:length(markerList)) {
  if(dim(marker[[i]])[1]>0){
    tmp=data.frame(matrix(ncol = 3,nrow = dim(marker[[i]])[1]))
    colnames(tmp)=c('cluster','gene','Log2FC')
    tmp$cluster=i
    tmp$gene=marker[[i]]$name
    tmp$Log2FC=marker[[i]]$Log2FC
    tmp=tmp[order(-tmp$Log2FC),]}
  if(dim(marker[[i]])[1]==0){
    tmp=data.frame(matrix(ncol = 3,nrow = 1))
    colnames(tmp)=c('cluster','gene','Log2FC')
    tmp$cluster=i
    tmp$gene='null'
    tmp$Log2FC='NA'}
  use<-rbind(use,tmp)}

write.csv(use,file = 'marker.csv',row.names = F)


