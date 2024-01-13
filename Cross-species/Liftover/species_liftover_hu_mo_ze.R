setwd('~/CH/CH-cross/Fig3/liftover/')
library(data.table)
#mm
mm9_peak_align=fread('./motif/mm9_align.bed')
mm10_peak=fread('./motif/mm10_peak.bed')

trans_mm9_mm10=data.frame('mm9'=paste(mm9_peak_align$V1,mm9_peak_align$V2,mm9_peak_align$V3,sep = '-'),'mm10'=paste(mm10_peak$V1,mm10_peak$V2,mm10_peak$V3,sep = '-'))
rownames(trans_mm9_mm10)=trans_mm9_mm10$mm10

mm10_peak_align=fread('./motif/mm10_align.bed')
mm10_peak_align=paste(mm10_peak_align$V1,mm10_peak_align$V2,mm10_peak_align$V3,sep = '-')
mm10_peak_align[1:5]
trans_mm9_mm10[mm10_peak_align,]$mm9->mm9_align

hg38_peak=fread('./motif/mm10toHg38_peak.bed')
hg38_peak=paste(hg38_peak$V1,hg38_peak$V2,hg38_peak$V3,sep = '-')
trans_mm10_hg38=data.frame(mm10=mm10_peak_align,hg38=hg38_peak)
#trans_mm10_hg38=trans_mm10_hg38[!duplicated(trans_mm10_hg38$hg38),]
#rownames(trans_mm10_hg38)=trans_mm10_hg38$hg38

mm9_peak_total=fread('./mm9_peak.bed')
head(mm9_peak_total)
mm9_peak_total$V4=paste(mm9_peak_total$V1,mm9_peak_total$V2,mm9_peak_total$V3,sep = '-')

mm9_peak_total$tag='unalign'
mm9_peak_total[mm9_peak_total$V4%in%mm9_align,]$tag='align'

table(mm9_peak_total$tag)

mm10_conservehuman=fread('./motif/mm10_conservehuman.bed')
mm10_conservehuman=paste(mm10_conservehuman$V1,mm10_conservehuman$V2,mm10_conservehuman$V3,sep = '-')
mm10_conservehuman=unique(mm10_conservehuman)

library(dplyr)
intersect(rownames(trans_mm10_hg38),mm10_conservehuman)%>%length()
mm10_conservehuman_1=trans_mm10_hg38[trans_mm10_hg38$hg38%in%mm10_conservehuman,]$mm10%>%unique()

intersect(mm10_conservehuman_1,rownames(trans_mm9_mm10))%>%length()
mm9_conservehuman=trans_mm9_mm10[mm10_conservehuman_1,]$mm9


mm9_peak_total[mm9_peak_total$V4%in%mm9_conservehuman,]$tag='conserve'
table(mm9_peak_total$tag)

setwd('~/CH/CH-cross/Fig3/liftover/')
dr11_peak_align=fread('./motif/dr11_align.bed')
dr11_peak_align=paste(dr11_peak_align$V1,dr11_peak_align$V2,dr11_peak_align$V3,sep = '-')
dr11_peak_align[1:5]


hg38_peak=fread('./motif/dr11toHg38_peak.bed')
hg38_peak=paste(hg38_peak$V1,hg38_peak$V2,hg38_peak$V3,sep = '-')
trans_dr11_hg38=data.frame(dr11=dr11_peak_align,hg38=hg38_peak)


dr11_peak_total=fread('./dr11_peak_0913.bed')
head(dr11_peak_total)
dr11_peak_total$V4=paste(dr11_peak_total$V1,dr11_peak_total$V2,dr11_peak_total$V3,sep = '-')

dr11_peak_total$tag='unalign'
dr11_peak_total[dr11_peak_total$V4%in%dr11_peak_align,]$tag='align'

table(dr11_peak_total$tag)

dr11_conservehuman=fread('./motif/dr11_conservehuman.bed')
dr11_conservehuman=paste(dr11_conservehuman$V1,dr11_conservehuman$V2,dr11_conservehuman$V3,sep = '-')
dr11_conservehuman=unique(dr11_conservehuman)

library(dplyr)

dr11_conservehuman_1=trans_dr11_hg38[trans_dr11_hg38$hg38%in%dr11_conservehuman,]$dr11%>%unique()


dr11_peak_total[dr11_peak_total$V4%in%dr11_conservehuman_1,]$tag='conserve'
table(dr11_peak_total$tag)

dat_all=readxl::read_xlsx('./liftover.xlsx')
head(dat_all)
colnames(dat_all)[1]='Species'
dat_all=melt(dat_all)

dat_all$Species=factor(dat_all$Species,levels = rev(c('human','mouse','zebrafish')))
dat_all$variable=factor(dat_all$variable,levels = c('unalign','align','conserve'))
dat_all$value=log10(dat_all$value+1)

p=ggplot(dat_all) +geom_bar(aes(x =Species , y = value, fill = variable), stat = "identity", 
           width = 0.75, color = "black") + theme_classic() +coord_flip()+
  theme(axis.text = element_text(size= 11, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(2, 4), 'lines'))+scale_fill_manual(values = c('#00AFBB', '#E7B800','#6495ED'))+ylab('Number of accessible regions(log10)')

p

ggsave(p,filename = '~/CH/CH-cross/Fig3/liftover/species_liftover_hu_mo_ze_new.pdf',height = 5,width = 10)
