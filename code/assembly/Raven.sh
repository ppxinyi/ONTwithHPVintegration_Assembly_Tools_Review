#!/bin/bash
#SBATCH --job-name=Raven_um47
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60gb
#SBATCH --time=100:00:00
#SBATCH --account=chadbren99
#SBATCH --partition=standard
#SBATCH --output=r.log
#SBATCH --error=r.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate raven-env

mkdir -p UM47
cd UM47
raven \
    --threads 10 \
    --polishing-rounds 0 \
    --graphical-fragment-assembly UM47_raven.gfa \
    /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/UM47/ONT-2-UM47_864-LP_1_BD/864-LP_1/20200917_2223_X4_FAO29991_fb3c45b4/fastq_pass/FAO29991_combined.fastq.gz \
    > UM47_raven.fasta


mkdir -p align_result
cd align_result
module load Bioinformatics
module load samtools
module load minimap2
ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Raven_um47/UM47/UM47_raven.fasta
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Flye/hg38p14_HPV16.fa
THREADS=14

minimap2 -d ${REF%.fa}.mmi $REF
minimap2 -t $THREADS -ax asm5 --cs=long -Y  ${REF%.fa}.mmi $ASM \
  | samtools sort -@ 8 -o UM47.bam
samtools index UM47.bam


BAM=UM47.bam
OUTDIR=work_hpv
TOOLDIR=/scratch/chadbren_root/chadbren99/ppxinyi/canu_example/artifical_seq/W_result/3level/Level1/work_hpv/minimap2-2.30_x64-linux
chmod +x $TOOLDIR/minimap2 $TOOLDIR/k8 $TOOLDIR/paftools.js
export PATH="$TOOLDIR:$PATH"
mkdir -p "$OUTDIR"

samtools view -h "$BAM" \
 | awk '$0 !~ /^@/ && $3=="HPV16"{print $1}' \
 | sort -u > "$OUTDIR/hpv.qnames.txt"

samtools view -@8 -N "$OUTDIR/hpv.qnames.txt" -b "$BAM" > "$OUTDIR/UM47.HPVreads.bam"
samtools index "$OUTDIR/UM47.HPVreads.bam"

samtools view -h "$OUTDIR/UM47.HPVreads.bam" \
  | $TOOLDIR/k8 $TOOLDIR/paftools.js sam2paf - \
  | sort -k1,1 -k3,3n > "$OUTDIR/UM47.HPVreads.paf"

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
"$OUTDIR/UM47.HPVreads.paf" > "$OUTDIR/UM47.hpv_chains.tsv"

