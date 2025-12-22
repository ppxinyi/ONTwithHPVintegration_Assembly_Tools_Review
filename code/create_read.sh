#!/bin/bash
#SBATCH --job-name=pbsim3no
#SBATCH --mail-user=ppxinyi@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=300gb
#SBATCH --time=100:00:00
#SBATCH --account=chadbren99
#SBATCH --partition=largemem
#SBATCH --output=1.log
#SBATCH --error=1.err   

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pbsim3

module load Bioinformatics
module load samtools
module load minimap2
module load bedtools2/2.31.1-zl7ag52

set -euo pipefail

# ---- inputs ----
INTEG_GZ="/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/L3/Level3.withHPV.fa.gz"
HOST="/nfs/turbo/oto-brenner-lab/Xinyi/MCPV/artifical_seq/prepare/hg38.fa"
BEDPE="/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/L3/100_1218planEFG.bedpe.tsv"
MODEL="/nfs/turbo/oto-brenner-lab/Xinyi/MCPV/artifical_seq/Example1/pbsim3/data/QSHMM-ONT.model"
SUFFIX="|HPVsiteLists"        
# W=30000             
LENM=10000; LENS=6000; LENMIN=1000; LENMAX=200000
INTEG="/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/L3/Level3.withHPV.fa"
# #---- prep ----
# INTEG="${INTEG_GZ%.gz}"; [[ -f $INTEG ]] || gunzip -c "$INTEG_GZ" > "$INTEG"
# samtools faidx "$INTEG"
# cut -f1 "${INTEG}.fai" > contigs_full.txt                                  
# awk -F'\t' '{name=$1; sub(/\|.*/,"",name); print name}' "${INTEG}.fai" > contigs_use.txt  

# # HOST subgroup
# samtools faidx "$HOST" -r contigs_use.txt -o only_used_contigs.fa


# # ---- breakpoint windows ----
# awk -v W=$W -v SUF="$SUFFIX" 'BEGIN{OFS="\t"} NR==1{next}
# {nm=($7==""?"JUNC"NR:$7); m=$2; s=m-W; if(s<0)s=0; e=m+int(W*1.4); print $1 SUF,s,e,nm"|A"}' \
#   "$BEDPE" > hpv_junc_windows.bed

# awk 'NR==FNR{ok[$1]=1;next} ok[$1]' contigs_full.txt hpv_junc_windows.bed > hpv_junc_windows.integ.bed
# bedtools getfasta -fi "$INTEG" -bed hpv_junc_windows.integ.bed -name -fo hpv_windows.fa

# # ---- PBSIM3 3 part ----
# mkdir -p out_integ out_host out_spike
# cd out_integ 
# pbsim --strategy wgs --method qshmm --qshmm "$MODEL" --genome "$INTEG" \
#   --depth 15 --length-mean $LENM --length-sd $LENS --length-min $LENMIN --length-max $LENMAX \
#   --seed 2101 
# cd ..
# cd out_host
# pbsim --strategy wgs --method qshmm --qshmm "$MODEL" --genome "$HOST" \
#   --depth 10  --length-mean $LENM --length-sd $LENS --length-min $LENMIN --length-max $LENMAX \
#   --seed 0901 
# cd ..
# cd out_spike
# pbsim --strategy wgs --method qshmm --qshmm "$MODEL" --genome /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/L3_reads/hpv_windows.fa \
#   --depth 5 --length-mean $LENM --length-sd $LENS --length-min $LENMIN --length-max $LENMAX \
#   --seed 3102 
# cd ..
conda activate seqkit_env
module load Bioinformatics
module load samtools
module load minimap2
module load bedtools2/2.31.1-zl7ag52
# ---- merge ----
# seqkit replace -p '^(\S+)' -r 'INTEG_${1}' out_integ/*.fq* > r_integ.fq
# seqkit replace -p '^(\S+)' -r 'HOST_${1}'  out_host/*.fq*  > r_host.fq
# # seqkit replace -p '^(\S+)' -r 'SPIKE_${1}' out_spike/*.fq* > r_spike.fq
# cat r_integ.fq r_host.fq | seqkit shuffle -s 42 | gzip > ONT.mix.fq.gz
# echo "[OK] ONT.mix.fq.gz"

# cat only_used_contigs.fa /nfs/turbo/oto-brenner-lab/Xinyi/MCPV/artifical_seq/HPV16.fa > host_plus_HPV16.fa
# minimap2 -d host_plus_HPV16.mmi host_plus_HPV16.fa
# samtools faidx host_plus_HPV16.fa
# awk '{print $1 "\t0\t" $2}' host_plus_HPV16.fa.fai > all.bed
# awk '/HPV/ {print $1 "\t0\t" $2}' host_plus_HPV16.fa.fai > hpv.bed

# minimap2 -t 14 -x map-ont -a -Y host_plus_HPV16.mmi ONT.mix.fq.gz \
#   | samtools sort -@8 -o ONT.mix.bam
# samtools index ONT.mix.bam

# # count：HPV sa + cplit SA
# samtools view -F 4 ONT.mix.bam | awk '$3~/HPV|HPV16|MCPV/{print $1}' | sort -u | wc -l > hpv_mapped.count.txt
# samtools view -h ONT.mix.bam \
# | awk 'BEGIN{hpv="HPV|HPV16|HPV18|MCPV"} $1~/^@/{next}
#        { if (match($0,/SA:Z:([^ \t]+)/,m)) {
#            isHPV=($3~hpv); isCHR=($3~/^chr/);
#            if ((isHPV && m[1]~/chr[0-9XYM]/) || (isCHR && m[1]~/(HPV|HPV16|HPV18|MCPV)/)) print $1
#          }}' | sort -u | wc -l > hpv_chr_split.count.txt

# ####N50 for reads with HPV
# samtools view -F 4 ONT.mix.bam | awk '$3 ~ /HPV/ {print $1}' > hpv_ids1.txt
# samtools view ONT.mix.bam | awk '$0 ~ /SA:Z:.*HPV/ {print $1}' > hpv_ids2.txt
# cat hpv_ids1.txt hpv_ids2.txt | sort -u > hpv_mapped.read_ids.txt
# samtools view -F 4 ONT.mix.bam \
# | awk 'FNR==NR{keep[$1]=1; next} ($1 in keep){ if(!( $1 in len )) len[$1]=length($10) } END{for(k in len) print len[k]}' \
#    hpv_mapped.read_ids.txt - \
# | sort -nr > hpv_read_lengths.txt

# awk '{L+=$1; a[NR]=$1} END{half=L/2; s=0; for(i=1;i<=NR;i++){s+=a[i]; if(s>=half){print a[i]; exit}}}' \
#    hpv_read_lengths.txt > hpv_mapped.N50.txt

#    # 1) extract HPV alignments read IDs（only primary alignment）
# samtools view -F 0x904 ONT.mix.bam \
#   | awk '$3 ~ /HPV/ {print $1}' | sort -u > hpv_read_ids.txt
# # note：0x904 = 0x800(补充，2048) + 0x100(次要，256) + 0x004(未比对，4)

# # 2) get hPV priamry alignment reads bam file
# samtools view -@10 -N hpv_read_ids.txt -b ONT.mix.bam > ONT.mix.HPVreads.bam
# samtools index ONT.mix.HPVreads.bam

# 3) “only HPV-reads 的mean depth”

# HPV-only mean depth per bp
samtools depth -a -b hpv.bed ONT.mix.bam \
| awk '{sum+=$3; n++} END{printf("HPV-only aligned mean depth = %.3fx\n", sum/n)}' \
> hpv_only.Pmean_depth.txt

HPV_CONTIG="HPV16"
HPV_LEN=7906

samtools view -F 0x900 ONT.mix.bam "$HPV_CONTIG" \
| awk -v L="$HPV_LEN" '
{
  cigar = $6
  aligned = 0
  while (match(cigar, /^([0-9]+)([MIDNSHP=X])/, m)) {
    n = m[1] + 0
    op = m[2]
    if (op ~ /[M=X]/) aligned += n   # 消耗 reference 的对齐列
    cigar = substr(cigar, RLENGTH + 1)
  }
  sum += aligned
}
END {
  if (L>0) printf("HPV effective coverage (CIGAR-based) = %.2fx\n", sum/L)
  else print "ERROR: HPV_LEN is 0"
}' > HPV_effective_coverage.txt


samtools depth -a -b hpv.bed ONT.mix.bam \
| awk '{if($3>0) covered++; total++}
       END{printf("HPV-only coverage = %.2f%%\n", 100*covered/total)}' \
> hpv_only.coverage.txt

############################################
# 4) integration windows：mean depth + coverage
############################################

WIN2=20000
HPV_TAG="HPV16"  

awk -v W="$WIN2" -v HPV="$HPV_TAG" 'BEGIN{OFS="\t"} NR==1{next}
{
  c1=$1; s1=$2+0; e1=$3+0;
  c2=$4; s2=$5+0;
  id = (NF>=7 && $7!="") ? $7 : ("JUNC"NR)

  if (c1 ~ HPV && c2 !~ HPV) {chrom=c2; m=s2}
  else if (c2 ~ HPV && c1 !~ HPV) {chrom=c1; m=s1}
  else {
    chrom=c1; m=s1
  }

  start=m-W; if(start<0) start=0
  end=m+W
  if(start>end){tmp=start; start=end; end=tmp}

  print chrom, start, end, id
}' "$BEDPE" > host_junc_windows.${WIN2}bp.bed

# 
cut -f1 host_plus_HPV16.fa.fai | sort -u > ref_contigs.txt
awk 'NR==FNR{ok[$1]=1;next} ok[$1]' ref_contigs.txt host_junc_windows.${WIN2}bp.bed \
  > host_junc_windows.${WIN2}bp.onref.bed

echo "[OK] windows BED: host_junc_windows.${WIN2}bp.onref.bed"

# --- (A) all reads：every windows coverage + mean depth ---
# bedtools coverage last 3 columns: bases_covered, length, fraction_covered
# mean and mean_depth
bedtools coverage -a host_junc_windows.${WIN2}bp.onref.bed -b ONT.mix.bam -mean \
| awk 'BEGIN{OFS="\t"}
       {
         # output：chrom start end id cov_frac mean_depth
         cov=$(NF-1); mean=$(NF);
         print $1,$2,$3,$4,cov,mean
       }' > windows.${WIN2}bp.depth_cov.ALL.tsv

echo "[OK] windows.${WIN2}bp.depth_cov.ALL.tsv"

# --- (B) only HPV-reads：every windows coverage + mean depth ---
bedtools coverage -a host_junc_windows.${WIN2}bp.onref.bed -b ONT.mix.HPVreads.bam -mean \
| awk 'BEGIN{OFS="\t"}
       {
         cov=$(NF-1); mean=$(NF);
         print $1,$2,$3,$4,cov,mean
       }' > windows.${WIN2}bp.depth_cov.HPVreads.tsv

echo "[OK] windows.${WIN2}bp.depth_cov.HPVreads.tsv"

# mean coverage / mean depth）
{
  echo -e "BAM\tmean_cov\tmean_depth"
  awk 'BEGIN{OFS="\t"} {c+=$5; d+=$6; n++} END{print "ONT.mix.bam", c/n, d/n}' windows.${WIN2}bp.depth_cov.ALL.tsv
  awk 'BEGIN{OFS="\t"} {c+=$5; d+=$6; n++} END{print "ONT.mix.HPVreads.bam", c/n, d/n}' windows.${WIN2}bp.depth_cov.HPVreads.tsv
} > windows.${WIN2}bp.depth_cov.summary.txt

echo "[OK] windows.${WIN2}bp.depth_cov.summary.txt"
