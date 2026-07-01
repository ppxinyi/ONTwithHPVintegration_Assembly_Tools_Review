#!/bin/bash
#SBATCH --job-name=shastalocal3
#SBATCH --mail-user=ppxinyi@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb
#SBATCH --time=100:00:00
#SBATCH --account=chadbren0
#SBATCH --partition=standard
#SBATCH --output=3.log
#SBATCH --error=3.err    

source ~/miniconda3/etc/profile.d/conda.sh


# awk 'NR%4==2 { gsub(/[Nn]/,"A"); } { print }' /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/pbsim3reads/localread/20000/l3/L3selected.fq > EFGreads.NtoA.fq
# zcat /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Ratio1to1/30depthHQ/EFGreads/ONT.mix.fq.gz  > EFG30x.fastq
conda activate shasta

shasta \
  --input  /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/assembly/local/shasta/EFGreads.NtoA.fq \
  --assemblyDirectory /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/assembly/local/shasta/20000/L3 \
  --threads 14 \
  --ReadGraph.strandSeparationMethod 1 \
  --config Nanopore-Sep2020 \
  --Reads.minReadLength 1000 \
  --MarkerGraph.minCoverage 1 \
  --MarkerGraph.minCoveragePerStrand 1 \
  --MarkerGraph.minEdgeCoverage 1 
#   --config Nanopore-R9.4.1-400bps \

mkdir -p align_L330
cd align_L330

module load Bioinformatics
module load samtools
module load minimap2
module load bedtools2/2.31.1-zl7ag52

ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/assembly/local/shasta/20000/L3/Assembly.fasta
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/1Review_paper/hg38_HPV16.fa
THREADS=14

minimap2 -d ${REF%.fa}.mmi $REF
minimap2 -t $THREADS -ax asm5 --cs=long -Y  ${REF%.fa}.mmi $ASM \
  | samtools sort -@ 8 -o L3.bam
samtools index L3.bam

BAM=L3.bam
OUTDIR=work_hpv
TOOLDIR=/scratch/chadbren_root/chadbren99/ppxinyi/canu_example/artifical_seq/W_result/3level/Level1/work_hpv/minimap2-2.30_x64-linux
chmod +x $TOOLDIR/minimap2 $TOOLDIR/k8 $TOOLDIR/paftools.js
export PATH="$TOOLDIR:$PATH"
mkdir -p "$OUTDIR"

samtools view -h "$BAM" \
 | awk '$0 !~ /^@/ && $3=="HPV16"{print $1}' \
 | sort -u > "$OUTDIR/hpv.qnames.txt"

samtools view -@8 -N "$OUTDIR/hpv.qnames.txt" -b "$BAM" > "$OUTDIR/L3.HPVreads.bam"
samtools index "$OUTDIR/L3.HPVreads.bam"

samtools view -h "$OUTDIR/L3.HPVreads.bam" \
  | $TOOLDIR/k8 $TOOLDIR/paftools.js sam2paf - \
  | sort -k1,1 -k3,3n > "$OUTDIR/L3.HPVreads.paf"

awk '
# ctg1    chr6:...(+ ) -> HPV16:...(-) -> chr6:...(+) ...
BEGIN{FS=OFS="\t"}
{
  q=$1; qS=$3; str=$5; t=$6; tS=$8; tE=$9;             
  seg = sprintf("%s:%d-%d(%s)", t, tS+1, tE, str);    
  if (q!=cur){
     if (chain!="") print cur, chain;
     cur=q; chain=seg;
  } else {
     chain = chain " -> " seg;
  }
}
END{ if (chain!="") print cur, chain }' \
"$OUTDIR/L3.HPVreads.paf" > "$OUTDIR/L3.hpv_chains.tsv"

