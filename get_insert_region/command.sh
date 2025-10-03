
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Chromosome sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

module load mysql 
module load Bioinformatics 
module load bedtools/v
# UCSC RepeatMasker table dump
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -e \
"select genoName, genoStart, genoEnd, repClass, repFamily, repName \
 from rmsk where genoName like 'chr%';" hg38 > hg38_rmsk.tsv

# Convert to BED
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' hg38_rmsk.tsv > hg38_rmsk.bed
grep "LINE" hg38_rmsk.bed > hg38_LINE.bed
grep "SINE" hg38_rmsk.bed > hg38_SINE.bed

# Make windows (100 kb example, adjust as needed)
bedtools makewindows -g hg38.chrom.sizes -w 100000 > hg38_100kb.bed

# Compute GC content per window
bedtools nuc -fi hg38.fa -bed hg38_100kb.bed > hg38_100kb_gc.tsv
awk 'BEGIN{OFS="\t"} NR>1 {if($5<0.35) print $1,$2,$3,"lowGC",$5; \
 else if($5>0.55) print $1,$2,$3,"highGC",$5}' hg38_10kb_gc.tsv > hg38_gc_classified.bed

mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -e \
"select chrom, chromStart, chromEnd, otherChrom, otherStart, otherEnd, fracMatch \
 from genomicSuperDups where chrom like 'chr%';" hg38 > hg38_segdup.tsv

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"SegDup",$7}' hg38_segdup.tsv > hg38_segdup.bed

# Overlap with LINE/SINE
bedtools intersect -a integration_sites.bed -b hg38_LINE.bed -wa -wb > integ_in_LINE.tsv
bedtools intersect -a integration_sites.bed -b hg38_SINE.bed -wa -wb > integ_in_SINE.tsv

# Overlap with GC bins
bedtools intersect -a integration_sites.bed -b hg38_gc_classified.bed -wa -wb > integ_in_GC.tsv

# Overlap with segmental duplications
bedtools intersect -a integration_sites.bed -b hg38_segdup.bed -wa -wb > integ_in_SegDup.tsv


########Find repeat segment
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A -e \
"select chrom, chromStart, chromEnd, otherChrom, otherStart, otherEnd, fracMatch \
from genomicSuperDups where chrom like 'chr%';" hg38 > hg38_segdup.tsv
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"SegDup",$7}' hg38_segdup.tsv > hg38_segdup.bed
bedtools intersect -a integration_sites.bed -b hg38_segdup.bed -wa -wb > sites_in_segdup.tsv
