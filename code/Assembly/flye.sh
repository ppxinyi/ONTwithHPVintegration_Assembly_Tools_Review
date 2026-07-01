#!/bin/bash
#SBATCH --job-name=F1-30
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=250gb
#SBATCH --time=100:00:00
#SBATCH --account=chadbren0
#SBATCH --partition=largemem
#SBATCH --output=fs.log
#SBATCH --error=fs.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate flye-env

flye \
  --nano-hq /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/pbsim3reads/level1/L1ONT.mix.fq.gz \
  --out-dir  l1assembly_result \
  --threads 14 \
  --genome-size 3.1g \
  --asm-coverage 200 \
  --min-overlap 3000 \


module load Bioinformatics
module load samtools
module load minimap2
ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/assembly/whole/flye/l1/l1assembly_result/assembly.fasta
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/1Review_paper/hg38_HPV16.fa
THREADS=14

mkdir -p align_results
cd align_results

minimap2 -d ${REF%.fa}.mmi $REF
minimap2 -t $THREADS -ax asm5 --cs=long -Y  ${REF%.fa}.mmi $ASM \
  | samtools sort -@ 8 -o l1.bam
samtools index l1.bam


BAM=l1.bam
OUTDIR=work_hpv
TOOLDIR=/scratch/chadbren_root/chadbren99/ppxinyi/canu_example/artifical_seq/W_result/3level/Level1/work_hpv/minimap2-2.30_x64-linux
chmod +x $TOOLDIR/minimap2 $TOOLDIR/k8 $TOOLDIR/paftools.js
export PATH="$TOOLDIR:$PATH"
mkdir -p "$OUTDIR"

samtools view -h "$BAM" \
 | awk '$0 !~ /^@/ && $3=="HPV16"{print $1}' \
 | sort -u > "$OUTDIR/hpv.qnames.txt"

samtools view -@8 -N "$OUTDIR/hpv.qnames.txt" -F 0x100 -b "$BAM" > "$OUTDIR/l1.HPVreads.bam"
samtools index "$OUTDIR/l1.HPVreads.bam"

###if include second align remove -f 0x100
samtools view -h "$OUTDIR/l1.HPVreads.bam" \
  | $TOOLDIR/k8 $TOOLDIR/paftools.js sam2paf - \
  | sort -k1,1 -k3,3n > "$OUTDIR/l1.HPVreads.paf"

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
"$OUTDIR/l1.HPVreads.paf" > "$OUTDIR/l1.hpv_chains.tsv"
