load('~/CH/CH-zf/datu/Fig2/markerGS.rda')

df=data.frame()
for (i  in names(markerList@listData)) {
  dapeak_num=dim(markerList@listData[[i]])[1]
  tmp=data.frame(ct=i,dapeak_num=dapeak_num)
  df=rbind(df,tmp)
}
write.csv(df,file = '~/CH/CH-zf/datu/Fig2/dapeak_num.csv',quote = F,row.names = F)

#jss score
num=c()
m.list=list()
for (i  in names(markerList@listData)) {
  if(!grepl('unknown',i)){
    dapeak_num=dim(markerList@listData[[i]])[1]
    num=c(num,dapeak_num)
    print(dapeak_num)
    m.list[[i]]=markerList@listData[[i]]
  }
}
median(num)
names(m.list)

idy.ls <- m.list

idy=list()
for(i in 1:length(idy.ls)){
  idy[[i]]= as.data.frame(matrix(ncol=7,nrow = length(idy.ls[[i]])))
  idy[[i]]=idy.ls[[i]]
}

all.idy <- data.frame(matrix(ncol=7,nrow = 0))
for (i in 1:length(idy)) {
  tmp <- as.data.frame(idy[[i]])
  all.idy <- rbind(all.idy,tmp)
}

all.idy <- as.data.frame(all.idy)


peak <- paste0(all.idy$seqnames,'_',all.idy$start,'_',all.idy$end)
peak[1:5]
peak1 <- peak[!duplicated(peak)]
peak1[1:5]


pmat1<-pmat[peak1,]


for (i in names(m.list)) {
  m.list[[i]]$peak=paste0(m.list[[i]]$seqnames,'_',m.list[[i]]$start,'_',m.list[[i]]$end)
  rownames(m.list[[i]])=m.list[[i]]$peak
}
daps=m.list



meta=proj2@cellColData
meta$barcodes=rownames(meta)
pmat_pertype = c()
all.anno=unique(meta$anno)

for(i in all.anno){ 
  selected.cells = meta[meta$anno == i,]$barcodes
  pmat.extracted = pmat1[,selected.cells]
  avg_expr = rowSums(pmat.extracted)/ncol(pmat.extracted)
  pmat_pertype = cbind(pmat_pertype, avg_expr)}


pmat_pertype = as.data.frame(pmat_pertype)
colnames(pmat_pertype) = all.anno
propmat=pmat_pertype
library(Matrix)

dp <- data.frame(matrix(nrow=0,ncol=2))
colnames(dp) <- c('V1','V2')
for (i in all.anno) {
  selected.cells = meta[meta$anno== i,]$barcodes
  pmat.extracted = pmat1[,selected.cells]
  avg_expr = median(colSums(pmat.extracted))
  dat=data.frame('V1'=i,'V2'=avg_expr)
  dp=rbind(dp,dat)
}

depthmat=dp

source('~/CH/code_practice/specificity_data/func.R')
print("Normalizing proportions...")
depthnorm = mean(depthmat[,2])/depthmat[,2]
logdepthnorm = mean(log10(depthmat[,2]))/log10(depthmat[,2])
propmatnormbylogdepth = t(t(propmat)*logdepthnorm)


print("Calculating specificity scores...")
marker_specificities_out = specificity_scorer(propmatnormbylogdepth)
markerdup = marker_specificities_out^2
markering = markerdup * propmatnormbylogdepth
rownames(markering) = rownames(propmat)
colnames(markering) = colnames(propmat)
markering[1:5,1:5]
markeringlong = melt(as.matrix(markering))
markeringlong = markeringlong[,c(3,1,2)]
colnames(markeringlong) = c("specificity_score","locusID","lineage")



#GO
#peak 2 gene

require(ChIPseeker)
require(clusterProfiler)
library(GenomicFeatures)
library(TxDb.Drerio.UCSC.danRer11.refGene)
getTxDbGenes <- function(txdb = NULL, orgdb = NULL, gr = NULL, ignore.strand = TRUE){
  
  if (is.null(genome)) {
    if (is.null(txdb) | is.null(orgdb)) {
      stop("If no provided genome then you need txdb and orgdb!")
    }
  }
  
  if (is.null(gr)) {
    genes <- GenomicFeatures::genes(txdb)
  }else {
    genes <- suppressWarnings(subsetByOverlaps(GenomicFeatures::genes(txdb), gr, ignore.strand = ignore.strand))
  }
  
  if (length(genes) > 1) {
    mcols(genes)$symbol <- suppressMessages(mapIds(orgdb, 
                                                   keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID", 
                                                   multiVals = "first"))
    genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
    names(genes) <- NULL
    out <- genes
  }else {
    out <- GRanges(seqnames(gr), ranges = IRanges(0, 0), gene_id = 0, symbol = "none")[-1]
  }
  
  return(out)
  
}
txdb <- TxDb.Drerio.UCSC.danRer11.refGene
library(org.Dr.eg.db)
genes <- getTxDbGenes(txdb=txdb,orgdb=org.Dr.eg.db)
id_symbol <- data.frame('geneId'=genes$gene_id,'symbol'=genes$symbol)



#GO heatmap
all.result <- as.data.frame(matrix(nrow=0,ncol = 3))
colnames(all.result)=c('Description','pvalue','lineage')
ct=names(m.list)
ct=c('Cardiomyocyte','Stromal cell','B cell','Cone photoreceptor cell','Rod photoreceptor cell','Granule cell','Olfactory neuron','Radial glia','Retina bipolar cell','Endothelial cell','Enterocyte','Intestinal bulb cell','Hepatocyte','Nephron cell','Epithelial cell','Ionocyte','Pancreatic cell')
#get bp result top5/all
for (i in ct) {
  message(i)
  lineage_re <- GRanges(m.list[[i]])
  peakAnno <- annotatePeak(lineage_re, tssRegion = c(-3000, 3000),
                           TxDb = txdb)
  peakAnno_df <- as.data.frame(peakAnno)
  peakAnno_df <- peakAnno_df[peakAnno_df$annotation!='Distal Intergenic',]
  
  peak_anno_df1 <- merge(peakAnno_df,id_symbol,by='geneId')
  peak_anno_df1<- peak_anno_df1[,c("seqnames", "start", "end", "symbol")]
  peak_anno_df1$peak=paste0(peak_anno_df1$seqnames,'_',peak_anno_df1$start,'_',peak_anno_df1$end)
  peak_anno_df1=peak_anno_df1[!duplicated(peak_anno_df1$peak),]
  
  peak_raw <- markeringlong[markeringlong$lineage==i,]
  peak_raw <- peak_raw[order(peak_raw$specificity_score,decreasing = T),]
  peak_raw$locusID=as.character(peak_raw$locusID)
  
  rownames(peak_anno_df1)=peak_anno_df1$peak
  peak_anno_df2=peak_anno_df1[peak_raw$locusID[1:500],]
  #peak_anno_df2=peak_anno_df1
  gene.use=peak_anno_df2$symbol
  gene.use <- gene.use[!is.na(gene.use)]
  gene.use <- gene.use[!duplicated(gene.use)]
  gene.use
  gene_list <- bitr(gene.use,fromType = "SYMBOL",toType =  "ENTREZID",
                    OrgDb =org.Dr.eg.db)
  go_list_bp<- enrichGO(gene = gene_list$ENTREZID,
                        OrgDb = org.Dr.eg.db, 
                        ont = 'BP',
                        pvalueCutoff =0.2, qvalueCutoff = 0.2,
                        readable = TRUE)
  
  go_list_bp@result[,c('Description','pvalue')]$Description[1:3]
  result <- go_list_bp@result[,c('Description','pvalue')]
  result$lineage=i
  result <- result[sort(result$pvalue,index.return=T)$ix[1:2],]
  all.result <- rbind(all.result,result)
}

top2 <- all.result

#heatmap 
top2$pvalue=-log10(top2$pvalue)
d <- top2%>% select("Description", "lineage", "pvalue") %>% dcast(Description ~ lineage)
d <- as.data.frame(d)
rownames(d)=d$Description
d=d[,-1]


bp.use <- unique(top3$Description)
d[is.na(d)]=0

d=d[bp.use,]
d=d[,ct]
mtx=d
library(pheatmap)
library(viridis)
p=pheatmap::pheatmap(mtx,
                     color = viridis(10,option = 'F'),
                     scale = 'none',
                     border_color = 'NA',
                     cluster_rows = F,
                     cluster_cols = F,
                     fontsize = 15,
                     width = 3,
                     height = 3)
p
ggsave(p,filename = '~/CH/CH-zf/datu/Fig2/Goterm/GO_heamap.pdf',height = 10,width = 20)