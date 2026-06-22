#!/bin/bash
#SBATCH --job-name=Raven_L330x
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=180gb
#SBATCH --time=100:00:00
#SBATCH --account=chadbren99
#SBATCH --partition=standard
#SBATCH --output=r.log
#SBATCH --error=r.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate raven-env

READS=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/pbsim3reads/L3ONT.mix.fq.gz


mkdir -p L330x
cd L330x
raven \
    --threads 14 \
    --polishing-rounds 1 \
    -u 4000 \
    -f 0.0005 \
    --graphical-fragment-assembly L330x_raven.gfa \
    $READS \
    > L330x_raven.fasta

mkdir -p align_result
cd align_result
module load Bioinformatics
module load samtools
module load minimap2

ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/assembly/whole/raven/L330x/L330x_raven.fasta
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/1Review_paper/hg38_HPV16.fa
THREADS=14

minimap2 -d ${REF%.fa}.mmi $REF
minimap2 -t $THREADS -ax asm5 --cs=long -Y  ${REF%.fa}.mmi $ASM \
  | samtools sort -@ 8 -o L330x.bam
samtools index L330x.bam



BAM=L330x.bam
OUTDIR=work_hpv
TOOLDIR=/scratch/chadbren_root/chadbren99/ppxinyi/canu_example/artifical_seq/W_result/3level/Level1/work_hpv/minimap2-2.30_x64-linux
chmod +x $TOOLDIR/minimap2 $TOOLDIR/k8 $TOOLDIR/paftools.js
export PATH="$TOOLDIR:$PATH"
mkdir -p "$OUTDIR"

samtools view -h "$BAM" \
 | awk '$0 !~ /^@/ && $3=="HPV16"{print $1}' \
 | sort -u > "$OUTDIR/hpv.qnames.txt"

samtools view -@8 -N "$OUTDIR/hpv.qnames.txt" -b "$BAM" > "$OUTDIR/L330x.HPVreads.bam"
samtools index "$OUTDIR/L330x.HPVreads.bam"

samtools view -h "$OUTDIR/L330x.HPVreads.bam" \
  | $TOOLDIR/k8 $TOOLDIR/paftools.js sam2paf - \
  | sort -k1,1 -k3,3n > "$OUTDIR/L330x.HPVreads.paf"

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
"$OUTDIR/L330x.HPVreads.paf" > "$OUTDIR/L330x.hpv_chains.tsv"

