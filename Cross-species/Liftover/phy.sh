#mm9 to mm10
~/CH/CH-cross/cluster/compare/liftOver mm9_align.bed ~/CH/CH-cross/cluster/liftover/mm9ToMm10.over.chain mm10_align.bed unmapped.bed
~/CH/CH-cross/cluster/compare/liftOver mm9_unalign.bed ~/CH/CH-cross/cluster/liftover/mm9ToMm10.over.chain mm10_unalign.bed unmapped.bed
~/CH/CH-cross/cluster/compare/liftOver mm9_conserve.bed ~/CH/CH-cross/cluster/liftover/mm9ToMm10.over.chain mm10_conserve.bed unmapped.bed


#phy
phy=~/CH/CH-cross/cluster/phylop/mm10.60way.phyloP60way.bw

# phyloP
for i in `ls mm10_*.bed`;do
tab_name=`echo $i | sed 's#.bed#.phytab#g'`
res_name=`echo $i | sed 's#.bed#.phyres#g'`
~/CH/CH-cross/cluster/compare/dr11/align/bigWigAverageOverBed $phy $i $tab_name 
awk '{print $1,$6}' $tab_name > $res_name
done

#dr11 to dr7
~/CH/CH-cross/cluster/compare/liftOver dr11_align.bed ~/CH/CH-cross/cluster/liftover/danRer11ToDanRer7.over.chain.gz dr7_align.bed unmapped.bed
~/CH/CH-cross/cluster/compare/liftOver dr11_unalign.bed ~/CH/CH-cross/cluster/liftover/danRer11ToDanRer7.over.chain.gz dr7_unalign.bed unmapped.bed
~/CH/CH-cross/cluster/compare/liftOver dr11_conserve.bed ~/CH/CH-cross/cluster/liftover/danRer11ToDanRer7.over.chain.gz dr7_conserve.bed unmapped.bed

~/CH/CH-cross/cluster/compare/liftOver dr11.bed ~/CH/CH-cross/cluster/liftover/danRer11ToDanRer7.over.chain.gz dr7.bed unmapped.bed
less dr11.bed|grep -vf unmapped.bed >dr11_align1.bed
phy=~/CH/CH-zf/datu/Fig2/phy/zf/vertebrate.phyloP8way.bw

# phyloP
for i in `ls dr7_*.bed`;do
tab_name=`echo $i | sed 's#.bed#.phytab#g'`
res_name=`echo $i | sed 's#.bed#.phyres#g'`
~/CH/CH-zf/datu/Fig2/phy/bigWigAverageOverBed $phy $i $tab_name 
awk '{print $1,$6}' $tab_name > $res_name
done