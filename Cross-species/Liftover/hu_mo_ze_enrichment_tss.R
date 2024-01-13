#tss enrichment
#mouse
library(readr)
setwd('~/CH/CH-cross/Fig3/liftover/motif/mouse_phy/')


phy_unalign <- read_delim('./mm10_unalign.phyres',delim = " ",col_names = F)
colnames(phy_unalign)=c('peak','phy_score')
phy_unalign$class='N'

phy_align <- read_delim('./mm10_align.phyres',delim = " ",col_names = F)
colnames(phy_align)=c('peak','phy_score')
phy_align$class='S'

phy_conserve <- read_delim('./mm10_conserve.phyres',delim = " ",col_names = F)
colnames(phy_conserve)=c('peak','phy_score')
phy_conserve$class='SF'

##unalign
peak=phy_unalign$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr
library(ArchR)
addArchRGenome('mm9')
tss=geneAnnoMm9$TSS

id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
unalign=length(id)/length(peak)

#align
peak=phy_align$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr
library(ArchR)
addArchRGenome('mm9')
tss=geneAnnoMm9$TSS

id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
align=length(id)/length(peak)

#conserve
peak=phy_conserve$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr
library(ArchR)
addArchRGenome('mm9')
tss=geneAnnoMm9$TSS

id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
conserve=length(id)/length(peak)
df=data.frame(class=c('N','S','SF'),enrich=c(unalign,align,conserve))
df$class=factor(df$class,levels = c('N','S','SF'))
m_df=df

#dr11
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(readr)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Drerio.UCSC.danRer11.refGene, OrgDb = org.Dr.eg.db)
tss=geneAnnotation$TSS
setwd('~/CH/CH-cross/Fig3/liftover/motif/zebrafish_phy/')
phy.ls=list.files(pattern = 'phyres')


phy_unalign <- read_delim('./dr7_unalign.phyres',delim = " ",col_names = F)
colnames(phy_unalign)=c('peak','phy_score')
phy_unalign$class='N'

phy_align <- read_delim('./dr7_align.phyres',delim = " ",col_names = F)
colnames(phy_align)=c('peak','phy_score')
phy_align$class='S'

phy_conserve <- read_delim('./dr7_conserve.phyres',delim = " ",col_names = F)
colnames(phy_conserve)=c('peak','phy_score')
phy_conserve$class='SF'


##unalign
peak=phy_unalign$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr


id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
unalign=length(id)/length(peak)

#align
peak=phy_align$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr

id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
align=length(id)/length(peak)

#conserve
peak=phy_conserve$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr


id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
conserve=length(id)/length(peak)
df=data.frame(class=c('N','S','SF'),enrich=c(unalign,align,conserve))
df$class=factor(df$class,levels = c('N','S','SF'))
zf_df=df

#human
library(readr)
setwd('~/CH/CH-cross/Fig3/liftover/motif/human_phy/')


phy_unalign <- read_delim('./hg38_unalign.phyres',delim = " ",col_names = F)
colnames(phy_unalign)=c('peak','phy_score')
phy_unalign$class='N'

phy_align <- read_delim('./hg38_align.phyres',delim = " ",col_names = F)
colnames(phy_align)=c('peak','phy_score')
phy_align$class='S'

phy_conserve <- read_delim('./hg38_conserve.phyres',delim = " ",col_names = F)
colnames(phy_conserve)=c('peak','phy_score')
phy_conserve$class='SF'

##unalign
peak=phy_unalign$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr
library(ArchR)
addArchRGenome('Hg38')
tss=geneAnnoHg38$TSS

id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
unalign=length(id)/length(peak)

#align
peak=phy_align$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr
library(ArchR)
addArchRGenome('Hg38')
tss=geneAnnoHg38$TSS

id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
align=length(id)/length(peak)

#conserve
peak=phy_conserve$peak
peak[1:5]
peak=unique(peak)
chr=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c1
start=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c2
end=reshape2::colsplit(peak,'-',names = c('c1','c2','c3'))$c3

peak.gr=GRanges(seqnames = chr,IRanges(start,end))
peak.gr
library(ArchR)
addArchRGenome('Hg38')
tss=geneAnnoHg38$TSS

id=queryHits(findOverlaps(peak.gr,tss))%>%unique()
conserve=length(id)/length(peak)
df=data.frame(class=c('N','S','SF'),enrich=c(unalign,align,conserve))
df$class=factor(df$class,levels = c('N','S','SF'))
h_df=df

h_df$species='Human'
zf_df$species='Zebrafish'
m_df$species='Mouse'
df=rbind(h_df,m_df,zf_df)
p=ggplot(df,aes(x =species, y = enrich,fill=class))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8,preserve = 'single'),
           width = 0.7,lwd=0.2) +scale_fill_manual(values = c('#00AFBB', '#E7B800','#6495ED'))+
  #coord_flip()+
  theme(legend.position = 'none')+ theme(axis.text = element_text(size= 11, colour = "black"),
                                         
                                         axis.title = element_text(size = 14, colour = "black", face = "bold"),                                                  panel.background = element_blank(),
                                         panel.grid = element_blank(),
                                         plot.margin = unit(rep(2, 4), 'lines'))+labs(y = "Enrichment relative of TSS", x = "")+theme_classic()+scale_y_continuous()
p
ggsave(p,filename = './hu_mo_ze_enrichment_tss.pdf',height = 15,width = 15)
