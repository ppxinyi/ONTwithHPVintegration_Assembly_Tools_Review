#!/bin/bash
#SBATCH --job-name=Next47Wgs
#SBATCH --account=chadbren0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=50g
#SBATCH --time=200:00:00
#SBATCH --output=2.out
#SBATCH --error=2.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ppxinyi@umich.edu


set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate nextdenovo_env

which python
python -c "import paralleltask; print('paralleltask ok')"


WORKDIR=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Nextdenover/UM47
cd "$WORKDIR"

/home/ppxinyi/miniconda3/envs/nextdenovo_env/bin/python /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Nextdenover/NextDenovo/nextDenovo run.cfg

cd Align

module load Bioinformatics
module load samtools
module load minimap2
module load bedtools2/2.31.1-zl7ag52

ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/0202shake/Assembly_results/Nextdenover/UM47/01_rundir/03.ctg_graph/nd.asm.fasta
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
