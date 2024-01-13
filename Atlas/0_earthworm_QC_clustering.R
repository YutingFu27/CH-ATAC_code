setwd('~/CH/CH-ew/bed_total/')


library(Seurat)
library(ArchR)
addArchRThreads(threads =18)
load('~/CH/CH-ew/res/ew_archR.rda')

library(data.table)


file.list=list.files(pattern='.bed.gz$')
sample.list=gsub('_ew_sort.bed.gz','',file.list)

ArrowFiles <- createArrowFiles(
  inputFiles =file.list,
  sampleNames=sample.list,
  minTSS =2,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  minFrags = 1000,
  genomeAnnotation = ew_GenomeAnnotation,
  geneAnnotation = ew_GeneAnnotation
)

ArrowFiles=list.files(pattern = '*.arrow')
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,geneAnnotation = ew_GeneAnnotation,genomeAnnotation = ew_GenomeAnnotation,copyArrows = F)
proj 


##rm doublets
doubScores <- addDoubletScores(
  input = proj1 ,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
proj1 <- filterDoublets(doubScores) 
proj1

df <- getCellColData(proj1, select = c("log10(nFrags)", "TSSEnrichment"))
ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment")+geom_hline(yintercept =2 ,lty='dashed',sie=0.25)+geom_vline(xintercept = log10(1000),lty='dashed',size=0.25)

##reduce dimension
res=1
addIterativeLSI_name=paste0('IterativeLSI_bin_res',res)
proj2 <- addIterativeLSI(
  ArchRProj = proj1,
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
saveRDS(proj2,file = 'proj2_bin.rds')

#peak-calling
proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = 'Clusters_bin_res1',force = TRUE)
pathToMacs2 <- '~/anaconda3/bin/macs2'
proj<- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = 'Clusters_bin_res1', 
  pathToMacs2=pathToMacs2,genomeSize=1.1e9,#dm genome size
)
proj2=addPeakMatrix(proj2)

#cell-peak 
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
pdf('peak.pdf',height = 10,width = 10)
plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = addClusters_name, embedding = addUMAP_name)
dev.off()

saveRDS(proj2,file = 'proj2_peak.rds')

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








