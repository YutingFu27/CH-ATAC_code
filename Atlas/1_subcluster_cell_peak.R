rm(list = ls())
library(data.table)
library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
#library(LSD)
library(reshape2)
addArchRThreads(threads = 6)
library(BSgenome.Dmelanogaster.UCSC.dm6)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Dmelanogaster.UCSC.dm6)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene, OrgDb = org.Dm.eg.db)


setwd('~CH/CH-dm/arrow/')

proj_total=readRDS('~CH/CH-dm/Fig1/proj2_peak_intergrate.rds')


lineage.ls=sort(unique(proj_total$lineage))
lineage.ls

for (i in lineage.ls) {
  message(i)
  #add info
  file_path=paste0('~CH/CH-dm/Fig1/1_subcluster/',i,'/')
  file=paste0(file_path,'proj2_bin.rds')
  proj2 <- readRDS(file)
  
  #peak calling
  addClusters_name='Clusters_bin_res1'
  
  
  proj2 <- addGroupCoverages(ArchRProj = proj2, groupBy = addClusters_name,force = TRUE)
  pathToMacs2 <- '/home/ggj/anaconda3/envs/SAMap/bin/macs2'
  proj2<- addReproduciblePeakSet(
    ArchRProj = proj2, 
    groupBy = addClusters_name, 
    pathToMacs2=pathToMacs2,
    genomeSize=1.2e8,#zf genome size,
    extendSummits = 25
  )
  proj2 <- addPeakMatrix(proj2)
  
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
  
  save_path=paste0('~CH/CH-dm/Fig1/1_subcluster/',i,'/')
  save_file=paste0(save_path,'proj2_peak.rds')
  saveRDS(proj2,file = save_file)
  
}


