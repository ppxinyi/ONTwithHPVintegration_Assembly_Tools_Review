#!/bin/bash
#SBATCH --job-name=w1
#SBATCH --mail-user=ppxinyi@umich.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --cpus-per-task=14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=180gb
#SBATCH --time=150:00:00
#SBATCH --account=chadbren99
#SBATCH --partition=standard
#SBATCH --output=s1.log
#SBATCH --error=s1.err
module load Bioinformatics
module load samtools
module load minimap2

mkdir -p l1
cd l1
READS=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/pbsim3reads/level1/L1ONT.mix.fq.gz
# # assemble long reads
/scratch/chadbren_root/chadbren99/ppxinyi/wtdbg2/wtdbg2 -x ont -g 2g -t 14 -i "$READS" -fo l1

# echo "step1 finish"
# # derive consensus
/scratch/chadbren_root/chadbren99/ppxinyi/wtdbg2/wtpoa-cns -t 14 -i l1.ctg.lay.gz -fo l1.fa

mkdir -p align_result
cd align_result
ASM=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/assembly/whole/wtdbg2/l1/l1.fa
REF=/scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/review_paper/1Review_paper/hg38_HPV16.fa
THREADS=14

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
 | awk '$3=="HPV16"{print $1}' \
 | sort -u > "$OUTDIR/hpv.qnames.txt"

samtools view -@8 -N "$OUTDIR/hpv.qnames.txt" -b "$BAM" > "$OUTDIR/l1.HPVreads.bam"
samtools index "$OUTDIR/l1.HPVreads.bam"

samtools view -h -F 0x100 "$OUTDIR/l1.HPVreads.bam" \
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


# conda activate seqkit_env
# # seqkit stats /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/Assembly_result/Flye/Level1/assembly_result/assembly.fasta > contigresult.txt
# seqkit grep -f /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/Assembly_result/Wbtgd2/Level1/align_result/work_hpv/hpv.qnames.txt /scratch/chadbren_root/chadbren99/ppxinyi/Artifical_seq_reviewpaper/12_18new/Assembly_result/Wbtgd2/Level1/l1.fa > l1.HPVcontigs.fa
# seqkit stats l1.HPVcontigs.fa > l1.HPVcontigs.stats.txt
