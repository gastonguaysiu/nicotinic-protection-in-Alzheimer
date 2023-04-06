#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=rrg-XXXX
#SBATCH --mem=64G

# APP nicotine bioinformatics pipeline

# Record the start time
date +"%T"

# Load necessary modules
module load sra-toolkit/2.8.2-1
module load bwa/0.7.17
module load samtools/1.9
module load bedtools/2.27.1
echo 'sra-toolkit, bwa, samtools, and bedtools have been loaded'

# Define path shortcuts
ps='/home/gsiu/scratch/scripts'
pg='/home/gsiu/scratch/genomes'
pp='/home/gsiu/scratch/projects'
px='/app_nic'

# Create an array with the fastq file codes
name=(
  'SRR4238596'
  'SRR4238597'
  'SRR4238598'
  'SRR4238599'
  'SRR4238600'
  'SRR4238601'
  'SRR4238602'
  'SRR4238603'
  'SRR4238604'
  'SRR4238605'
  'SRR4238606'
  'SRR4238607'
  'SRR4238608'
)

# Download raw data
echo 'getting the raw data'

# Create a directory for raw files if it does not exist
raw_dir="$pp$px/raw"
if [ ! -d "$raw_dir" ]; then
  mkdir "$raw_dir"
fi

# Download fastq files for each item in the name array using sra-toolkit
for i in "${name[@]}"; do
  fastq-dump --split-3 --clip --dumpbase --skip-technical -O $pp$px/raw/ $i
  echo "$i fastq downloaded"
  echo ""
done

date +"%T"

# Index bwa for the mouse reference genome (mm10)
bwa index $pg/mm10_ref/mm10_ref.fa
date +"%T"

# Create sam files
echo 'on to sam files'

# Create a directory for sam files if it does not exist
sam_dir="$pp$px/sam"
if [ ! -d "$sam_dir" ]; then
  mkdir "$sam_dir"
fi

# Convert fastq to sam for each item in the name array
for i in "${name[@]}"; do
  pe1="_1.fastq"
  pe2="_2.fastq"
  temp1="$i$pe1"
  temp2="$i$pe2"

  echo ""
  echo "starting $i sam"
  echo ""
  bwa mem $pg/mm10_ref/mm10_ref.fa $pp$px/raw/$temp1 $pp$px/raw/$temp2 > $pp$px/sam/$i.sam
  echo ""
  echo "$i sam complete"
  echo ""
done

date +"%T"

# Create bam files
echo 'lets make bam files'

# Create a directory for bam files if it does not exist
bam_dir="$pp$px/bam"
if [ ! -d "bam_dir" ]; then
  mkdir $pp$px/bam
fi

# Convert sam to bam for each item in the name array
for i in "${name[@]}"; do
  samtools view -Sb $pp$px/sam/$i.sam > $pp$px/bam/$i.bam
  echo ""
  echo "$i bam complete"
  echo ""
done

date +"%T"

# Create bed files
echo 'the final step is bed files'

# Create a directory for bed files if it does not exist
bed_dir="$pp$px/bed"
if [ ! -d "$bed_dir" ]; then
  mkdir "$bed_dir"
fi

# Convert bam to bed for each item in the name array
for i in "${name[@]}"; do
  bedtools bamtobed -i $pp$px/bam/$i.bam > $pp$px/bed/$i.bed
  echo ""
  echo "$i bed complete"
  echo ""
done

date +"%T"

# Sort sam files
echo 'on to sorting of sam files'

# Create a directory for sorted sam files if it does not exist
s_sam_dir="$pp$px/s_sam"
if [ ! -d "$s_sam_dir" ]; then
  mkdir "$s_sam_dir"
fi

# Sort sam files for each item in the name array
for i in "${name[@]}"; do
  echo ""
  echo "starting to sort $i sam"
  echo ""
  sort -k 3,3 -k 4,4n $pp$px/sam/$i.sam > $pp$px/s_sam/$i.sam
  echo ""
  echo "$i sam is now sorted"
  echo ""
done

date +"%T"

# Process sorted sam files
echo 'processing sorted sam files...'

# Create a directory for processed files if it does not exist
process1_dir="$pp$px/process1"
if [ ! -d "$process1_dir" ]; then
  mkdir "$process1_dir"
fi

# Assemble expressed genes and transcripts for each item in the name array
for i in "${name[@]}"; do

  # Create a directory for each output if it does not exist
  if [ ! -d "$pp$px/process1/$i" ]; then
    mkdir "$pp$px/process1/$i"
  fi

  # Run cufflinks to process the sorted sam file and update assemblies.txt
  echo ""
  echo "starting to process $i"
  echo ""
  $pt/cufflinks/cufflinks -o $pp$px/process1/$i $pp$px/s_sam/$i.sam
  echo ""
  echo "$i sorted sam is now processed"
  echo "adding $i path and name to assemblies.txt"
  echo "$pp$px/process1/$i/transcripts.gtf" >> $pp$px/assemblies.txt
done

date +"%T"

# Merge assemblies using cuffmerge
echo 'merging assemblies...'

# Create a directory for merged assemblies if it does not exist
merged_asm_dir="$pp$px/merged_asm"
if [ ! -d "$merged_asm_dir" ]; then
  mkdir "$merged_asm_dir"
fi

# Run cuffmerge to create a single merged transcriptome annotation
echo ""
echo "creating a single merged transcriptome annotation"
echo ""
$pt/cufflinks/cuffmerge -g $pg/gtf_gff/mouse/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf -o $pp$px/merged_asm $pp$px/assemblies.txt
echo ""
echo "complete merged transcriptome annotation"
echo ""
date +"%T"

# Quantify gene and transcript expression in RNA-Seq samples
echo 'quantifying transcript expression...'

# Create a directory for quantification results if it does not exist
quant_cxb_dir="$pp$px/quant_cxb"
if [ ! -d "$quant_cxb_dir" ]; then
  mkdir "$quant_cxb_dir"
fi

# Run
# Quantify gene and transcript expression in RNA-Seq samples
echo 'quantifying transcript expression...'

# Create a directory for quantification results if it does not exist
quant_cxb_dir="$pp$px/quant_cxb"
if [ ! -d "$quant_cxb_dir" ]; then
  mkdir "$quant_cxb_dir"
fi

# Run cuffquant for each item in the name array
for i in "${name[@]}"; do
  echo "working in $i"
  $pt/cufflinks/cuffquant -o $pp$px/quant_cxb -b $pg/mm10_ref/mm10ref.fa $pg/gtf_gff/mouse/gencode.vM23.chr_patch_hapl_scaff.annotation.gtf $pp$px/s_sam/$i.sam
  echo ""
  echo "$i is now quantified in a cxb file"
  echo ""
done

date +"%T"

