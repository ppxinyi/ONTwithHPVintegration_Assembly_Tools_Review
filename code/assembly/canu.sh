#!/bin/bash
#SBATCH --job-name=canu47WGS
#SBATCH --mail-user=ppxinyi@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=180gb
#SBATCH --time=200:00:00
#SBATCH --account=chadbren0
#SBATCH --partition=standard
#SBATCH --output=ca1.log
#SBATCH --error=ca1.err

module load gcc/11.2.0
module load Bioinformatics
module load samtools
export PATH=/scratch/chadbren_root/chadbren99/ppxinyi/canu-2.3/bin:$PATH
canu -version
which canu

# canu -p human_mcpvCB2 -d CB2TargetCanu \
#   genomeSize=0.5m \
#   -nanopore-raw  '/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Flye/UM47_local_reads/selected.fq' \
#   useGrid=false maxThreads=12 maxMemory=100 \
#   minReadLength=1000 \
#   corOutCoverage=60 \
#   corMhapSensitivity=high corMaxEvidenceErate=0.30 \
#   correctedErrorRate=0.105 \
#   merylMemory=34g  merylThreads=12 \
#   stopOnLowCoverage=2 minInputCoverage=2 \
#   ovsMemory=26g \
#   batMemory=90g  batThreads=12 \
#   corMhapThreads=3 corMhapConcurrency=4 \
#   saveReads=true bamOutput=false

/scratch/chadbren_root/chadbren99/ppxinyi/canu-2.3/bin/canu -p UM47 -d UM47 \
  genomeSize=1m \
  -nanopore /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/UM47/ONT-2-UM47_864-LP_1_BD/864-LP_1/20200917_2223_X4_FAO29991_fb3c45b4/fastq_pass/FAO29991_combined.fastq.gz \
  useGrid=false maxThreads=14 maxMemory=100 \
  minReadLength=2000 \
  corOutCoverage=60 \
  corMhapSensitivity=normal corMaxEvidenceErate=0.30 \
  correctedErrorRate=0.120 \
  merylMemory=34g  merylThreads=14 \
  ovsMemory=36g \
  batMemory=60g  batThreads=14 \
  saveReads=true bamOutput=false \
  stopOnLowCoverage=2 minInputCoverage=2


module load Bioinformatics
module load samtools
module load minimap2
module load bedtools2/2.31.1-zl7ag52

ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Canu/UM47/UM47.contigs.fasta
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Flye/hg38p14_HPV16.fa
THREADS=14

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
