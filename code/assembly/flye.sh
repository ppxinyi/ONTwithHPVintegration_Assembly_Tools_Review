#!/bin/bash
#SBATCH --job-name=F2
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=400gb
#SBATCH --time=100:00:00
#SBATCH --account=chadbren99
#SBATCH --partition=largemem
#SBATCH --output=fs222.log
#SBATCH --error=fs222.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate flye-env

flye \
  --nano-raw /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/L1_reads/ONT.mix.fq.gz \
  --out-dir  assembly_result \
  --threads 18 \
  --genome-size 2g \
  --asm-coverage 200 \
  --min-overlap 1000 \
  --scaffold 
mkdir -p align_result
cd align_result

module load Bioinformatics
module load samtools
module load minimap2
ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/Assembly_result/Flye/Level1/assembly_result/assembly.fasta
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/L1_reads/host_plus_HPV16.fa
THREADS=18

minimap2 -d ${REF%.fa}.mmi $REF
minimap2 -t $THREADS -ax asm5 --cs=long -Y  ${REF%.fa}.mmi $ASM \
  | samtools sort -@ 8 -o L1.bam
samtools index L1.bam


BAM=L1.bam
OUTDIR=work_hpv
TOOLDIR=/scratch/chadbren_root/chadbren99/ppxinyi/canu_example/artifical_seq/W_result/3level/Level1/work_hpv/minimap2-2.30_x64-linux
chmod +x $TOOLDIR/minimap2 $TOOLDIR/k8 $TOOLDIR/paftools.js
export PATH="$TOOLDIR:$PATH"
mkdir -p "$OUTDIR"

samtools view -h "$BAM" \
 | awk '$0 !~ /^@/ && $3=="HPV16"{print $1}' \
 | sort -u > "$OUTDIR/hpv.qnames.txt"

samtools view -@8 -N "$OUTDIR/hpv.qnames.txt" -b "$BAM" > "$OUTDIR/L1.HPVreads.bam"
samtools index "$OUTDIR/L1.HPVreads.bam"

samtools view -h "$OUTDIR/L1.HPVreads.bam" \
  | $TOOLDIR/k8 $TOOLDIR/paftools.js sam2paf - \
  | sort -k1,1 -k3,3n > "$OUTDIR/L1.HPVreads.paf"

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
"$OUTDIR/L1.HPVreads.paf" > "$OUTDIR/L1.hpv_chains.tsv"
