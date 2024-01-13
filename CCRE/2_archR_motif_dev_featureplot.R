rm(list = ls())
setwd("~/CH/CH-zf/datu/bed/")


library(Seurat)
library(ArchR)
addArchRChrPrefix(chrPrefix = FALSE)
library(reshape2)
addArchRThreads(threads =12)
library(BSgenome.Drerio.UCSC.danRer11)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Drerio.UCSC.danRer11)
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Drerio.UCSC.danRer11.refGene, OrgDb = org.Dr.eg.db)

proj2 <- readRDS('~/CH/CH-zf/datu/Fig2/JSD/proj2_coacc.rds')

proj2

#add motif anno
proj2 <- addMotifAnnotations(ArchRProj = proj2, motifSet = "encode", name = "Motif",force = TRUE)
load('~/CH/CH-zf/datu/Fig2/markerGS.rda')

#enrich motif
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersGS,
  ArchRProj = proj2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)


#add motif matrix to arrow
proj2 <- addBgdPeaks(proj2)

proj2 <- addDeviationsMatrix(
  ArchRProj = proj2, 
  peakAnnotation = "Motif",
  force = TRUE
)
saveRDS(proj2,file = '~/CH/CH-zf/datu/Fig2/homer/proj2_motif.rds')


proj2 <- readRDS('~/CH/CH-zf/datu/Fig2/homer/proj2_motif.rds')

##Motif deviatition
plotVarDev <- getVarDeviations(proj2, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 10, height = 10, ArchRProj = proj2, addDOC = FALSE)


motifs=c('HNF4a.NR','PU.1.ETS','CRX.Homeobox','RFX.HTH')
markerMotifs <- getFeatures(proj2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs=markerMotifs[c(3,2,4,1)]

#proj2@reducedDims$IterativeLSI_peak_res1
#proj2=addImputeWeights(proj2,reducedDims = 'IterativeLSI_peak_res1')
table(proj2$anno)

proj3=proj2[proj2$anno%in%c('Enterocyte','Intestinal bulb cell','Hepatocyte','B cell','T cell','Macrophage','Cone photoreceptor cell','Rod photoreceptor cell','GABAergic neuron','Oligodendrocyte','Primordial germ cell','Spermatocyte')]
proj3=addImputeWeights(proj3,reducedDims = 'IterativeLSI_peak_res1')

p <- plotGroups(ArchRProj = proj3, 
                groupBy = "anno", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj3)
)


ggsave(p[[1]],filename = '~/CH/CH-zf/datu/Fig2/motif/motif_dev1.pdf',height = 15,width = 15)
ggsave(p[[2]],filename = '~/CH/CH-zf/datu/Fig2/motif/motif_dev2.pdf',height = 15,width = 15)
ggsave(p[[3]],filename = '~/CH/CH-zf/datu/Fig2/motif/motif_dev3.pdf',height = 15,width = 15)
ggsave(p[[4]],filename = '~/CH/CH-zf/datu/Fig2/motif/motif_dev4.pdf',height = 15,width = 15)


#motif featureplot
motif <- c("PU.1", "HNF4A", "Ets1")
markerATAC <- getFeatures(proj2, select = paste(motif, collapse="|"), useMatrix = "MotifMatrix")
markerATAC <- sort(grep("z:", markerATAC, value = TRUE))
markerATAC=markerATAC[c(2,4,5)]
proj2=addImputeWeights(proj2,reducedDims = 'IterativeLSI_peak_res1')
p <- plotEmbedding(
  ArchRProj = proj2, 
  colorBy = "MotifMatrix", 
  name = markerATAC, 
  embedding = "UMAP_peak_res1",
  imputeWeights = getImputeWeights(proj2)
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
