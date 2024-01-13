#subcluster in cell-bin for drosophila
rm(list = ls())
library(data.table)
library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
#library(LSD)
library(reshape2)
addArchRThreads(threads = 8)
library(BSgenome.Dmelanogaster.UCSC.dm6)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Dmelanogaster.UCSC.dm6)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Dmelanogaster.UCSC.dm6.ensGene, OrgDb = org.Dm.eg.db)


setwd('~/CH/CH-dm/arrow/')

proj2=readRDS('~/CH/CH-dm/Fig1/proj2_peak_intergrate.rds')
table(proj2$lineage)
meta_total=as.data.frame(proj2@cellColData)

lineage.ls=sort(unique(proj2$lineage))
lineage.ls

for (i in lineage.ls) {
  message(i)
  #add info
  file_path=paste0('~/CH/CH-dm/Fig1/1_subcluster/',i,'/')
  file=paste0(file_path,'proj.rds')
  proj <- readRDS(file)
  meta=meta_total[proj$cellNames,]
  proj$linegae=i
  proj$subanno=meta$subanno
  
  #cell-bin cluster
  res=1
  addIterativeLSI_name=paste0('IterativeLSI_bin_res',res)
  proj2 <- addIterativeLSI(
    ArchRProj = proj,
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
  ##clustering
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
  
  #save plot&RDS
  plot_path=paste0('~/CH/CH-dm/Fig1/1_subcluster/',i,'/')
  plot=paste0(plot_path,'bin.pdf')
  pdf(plot,height = 10,width = 10)
  plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = addClusters_name, embedding = addUMAP_name)
  dev.off()
  
  save_path=paste0('~/CH/CH-dm/Fig1/1_subcluster/',i,'/')
  save_file=paste0(save_path,'proj2_bin.rds')
  saveRDS(proj2,file = save_file)
}
