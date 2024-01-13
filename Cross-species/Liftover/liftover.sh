#dr11 to hg38
~/CH-cross/cluster/compare/liftOver -minMatch=0.1 ../dr11_peak_0913.bed ~/CH-cross/cluster/liftover/danRer11ToHg38.over.chain.gz dr11toHg38_peak.bed unmapped.bed
less ../dr11_peak_0913.bed |grep -vf unmapped.bed > dr11_align.bed


bedtools intersect -f 0.5 -a ./dr11toHg38_peak.bed  -b ../hg38_peak_0918.bed -wa >dr11_conservehuman.bed


#hg38 to mouse
~/CH-cross/cluster/compare/liftOver -minMatch=0.1 ../hg38_peak_0918.bed ~/CH-cross/cluster/liftover/hg38ToMm10.over.chain.gz hg38Tomm10_peak.bed unmapped.bed
less ../hg38_peak_0918.bed |grep -vf unmapped.bed > hg38_align.bed
bedtools intersect -f 0.5 -a hg38Tomm10_peak.bed  -b mm10_peak.bed -wa >hg38_conservemouse.bed