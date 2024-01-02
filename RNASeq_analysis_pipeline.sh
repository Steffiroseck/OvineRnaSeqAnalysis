# Software Requirements
# FastQC v0.12.1 
# cutadapt version 1.15
# TrimmomaticPE
# HISAT2 version 2.1.0
# STAR
# featureCounts v2.0.1
# StringTie 2.2.1

# Path to the working directory where fastq files are present
cwd = "/home/afbi-roses@nigov.net/Sheep/"
cd $cwd

# Create directories to store outputs
mkdir 1.fatsqc
mkdir 2.trimmomatic
mkdir 3.hisat2
mkdir 3.star
mkdir 4.featurecounts
mkdir 4.stringtie
mkdir logs

# Step 1: Quality Control
# Run FASTQC (8 minutes per sample, 4 minutes per read)
mkdir 1.Fastqc.out
for R1 in 24.samples/*.fastq; 
do 
name=$(echo $R1|sed 's/.fastq//'|sed 's/24.samples\///'); 
echo "fastqc -o 1.Fastqc.out $R1 &> logs/$name.fastqclog.txt";
done > 1.fastqc.commands

# Step 2: Adapter trimming
# Run Trimmomatic
threads="8" #higher the number faster the speed; but make sure you have enough CPU support

for R1 in 24_Samples/*_R1_001.fastq.gz ; do
  R2="${R1%_R1_001.fastq.gz}_R2_002.fastq.gz"
  sample=$(echo $R1|sed 's/_R1_001.fastq.gz//'|sed 's/24_Samples\///'); 
  echo TrimmomaticPE -threads $threads "$R1" "$R2" 2. trimmomatic/"${sample}_R1_paired.fastq.gz" 2. trimmomatic/"${sample}_R1_unpaired.fastq.gz" 2. \
  trimmomatic/"${sample}_R2_paired.fastq.gz" 2. trimmomatic/"${sample}_R2_unpaired.fastq.gz" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 
done

# Step 3: Alignment to the reference genome using HISAT2
# 3.1: Build reference genome index. The genome should be in unzipped format or else it wont work
hisat2-build  Reference_geneome/ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna Reference_geneome/ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna

# 3.2 Now Align to the reference genome
for R1 in 2.trimmomatic/*_R1_paired.fastq.gz; do
R2=$(echo $R1| sed 's/_R1_/_R2_/'); 
sample=$(echo $R1|sed 's/_R1_paired.fastq.gz//'|sed 's/2.trimmomatic\///'); 
echo hisat2 -q --time --novel-splicesite-outfile 3.hisat2/$sample.tsv --summary-file 3.hisat2/$sample.txt \
--met-file 3.hisat2/$sample.txt --threads $threads -x Reference_geneome/ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna \
-1 $R1 -2 $R2 | tee >(samtools flagstat - > 3.hisat2/$sample.flagstat) \
| samtools sort -O BAM | tee 3.hisat2/$sample.bam \
| samtools index - 3.hisat2/$sample.bam.bai &> logs/$sample.hisat2.txt;
done

# Step 4: Alignment using STAR
# 4.1 Build genome index files
fastafile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna"
gtffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gtf"
threads=8

STAR --runThreadN $threads --runMode genomeGenerate --genomeDir StarIndex/ --genomeFastaFiles $fastafile --sjdbGTFfile $gtffile --sjdbOverhang 100

# 4.2 After indexing, we go for STAR alignment.The input files of STAR can be single-end or pair-end fastq files. If the annotation file (.gtf file) is provided, the accuracy of alignment can be increased.
indexDir='/mnt/sda1/RNA/40-815970407/Sheep/StarIndex/'
threads=8

for prefix in $(ls 2.trimmomatic/*.trim.fastq.gz | rev | cut -c 18- | rev | uniq)
do
 STAR --runThreadN $threads --genomeDir $indexDir  --readFilesIn "${prefix}.R1.trim.fastq.gz" "${prefix}.R2.trim.fastq.gz" \
 --outFileNamePrefix ${prefix}  --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --limitBAMsortRAM 1756793772
done

# Step 5 Generate read counts matrix using featureCounts
# When you want to analyze the data for differential gene expression analysis, it would be convenient to have counts for all samples in a single file (gene count matrix).
gtffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gtf"
featureCounts -T 8 -t 'gene' -g 'gene_id' -f -a $gtffile -o 4.featurecounts/LambAllSamples.featureCounts 3.hisat2/*.bam

## Since the featureCounts output has additional columns with information about genomic coordinates, gene length etc., 
## we can use the cut command to select only those columns that you are interested in. Columns 1 and sample wise counts columns
cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32 4.featurecounts/LambAllSamples.featureCounts > 4.featurecounts/LambAllSamples.featureCounts.Rmatrix

# STep 6 Generate read count matrix using stringtie
