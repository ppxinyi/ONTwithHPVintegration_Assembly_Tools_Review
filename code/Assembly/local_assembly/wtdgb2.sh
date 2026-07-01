#!/bin/bash
#SBATCH --job-name=w2-2
#SBATCH --mail-user=ppxinyi@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60gb
#SBATCH --time=150:00:00
#SBATCH --account=remills0
#SBATCH --partition=standard
#SBATCH --output=s2.log
#SBATCH --error=s2.err
module load Bioinformatics
module load samtools
module load minimap2

mkdir -p L2_test
cd L2_test

READS=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/pbsim3reads/localread/15000/l2/L2selected.fq
# # assemble long reads
/scratch/chadbren_root/chadbren99/ppxinyi/wtdbg2/wtdbg2 -x ont -X 100 -g 10m -e 2 -L 1500 -t 14 -i "$READS" -fo L2

# echo "step1 finish"
# # derive consensus
/scratch/chadbren_root/chadbren99/ppxinyi/wtdbg2/wtpoa-cns -t 14 -i L2.ctg.lay.gz -fo L2.fa

mkdir -p align_result
cd align_result
ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/assembly/local/Wdtgb2/15000/L2_test/L2.fa
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/1Review_paper/hg38_HPV16.fa
THREADS=14

minimap2 -d ${REF%.fa}.mmi $REF
minimap2 -t $THREADS -ax asm5 --cs=long -Y  ${REF%.fa}.mmi $ASM \
  | samtools sort -@ 8 -o L2.bam
samtools index L2.bam

BAM=L2.bam
OUTDIR=work_hpv
TOOLDIR=/scratch/chadbren_root/chadbren99/ppxinyi/canu_example/artifical_seq/W_result/3level/Level1/work_hpv/minimap2-2.30_x64-linux
chmod +x $TOOLDIR/minimap2 $TOOLDIR/k8 $TOOLDIR/paftools.js
export PATH="$TOOLDIR:$PATH"
mkdir -p "$OUTDIR"

samtools view -h "$BAM" \
 | awk '$3=="HPV16"{print $1}' \
 | sort -u > "$OUTDIR/hpv.qnames.txt"

samtools view -@8 -N "$OUTDIR/hpv.qnames.txt" -b "$BAM" > "$OUTDIR/L2.HPVreads.bam"
samtools index "$OUTDIR/L2.HPVreads.bam"

samtools view -h -F 0x100 "$OUTDIR/L2.HPVreads.bam" \
  | $TOOLDIR/k8 $TOOLDIR/paftools.js sam2paf - \
  | sort -k1,1 -k3,3n > "$OUTDIR/L2.HPVreads.paf"

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
"$OUTDIR/L2.HPVreads.paf" > "$OUTDIR/L2.hpv_chains.tsv"


# conda activate seqkit_env
# # seqkit stats /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/Assembly_result/Flye/Level1/assembly_result/assembly.fasta > contigresult.txt
# seqkit grep -f /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/Assembly_result/Wbtgd2/Level1/align_result/work_hpv/hpv.qnames.txt /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/Assembly_result/Wbtgd2/Level1/L2.fa > L2.HPVcontigs.fa
# seqkit stats L2.HPVcontigs.fa > L2.HPVcontigs.stats.txt