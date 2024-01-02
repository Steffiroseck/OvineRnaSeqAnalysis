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

# Step 1: Quality Control
# Run FASTQC (8 minutes per sample, 4 minutes per read)
mkdir 1.Fastqc.out
for R1 in 24.samples/*.fastq; 
do 
name=$(echo $R1|sed 's/.fastq//'|sed 's/24.samples\///'); 
echo "fastqc -o 1.Fastqc.out $R1 &> Logs/$name.fastqclog.txt";
done > 1.fastqc.commands

# Step 2: Adapter trimming
# Run Trimmomatic
#For paired-end data. Give the sample names in fastqfiler1 and fastq filer 2.
trimmomatic_path="/usr/bin/TrimmomaticPE"
threads="6" #higher the number faster the speed; but make sure you have enough CPU support

#Adapter clipping:
#Matches to potential adapter sequences at the end of reads, then removes reads shorter than 20 bases. 

#For the next step, copy the line as many times as the number of samples, then 
#add the sequences of your adapters after the -a flag (e.g. -a ATTGGCTTTGGGCAT), as well as changing the sample names:
#The sequences of Illumina's TruSeq adapters are proprietary information, but can be recieved by emailing Illumina customer support at info@illumina.com 

for R1 in *_R1_001.fastq.gz*
do
   R2=${R1//_R1_001.fastq.gz/_R2_001.fastq.gz}
   R1paired=${R1//.fastq.gz/_paired.fastq.gz}
   R1unpaired=${R1//.fastq.gz/_unpaired.fastq.gz}	
   R2paired=${R2//.fastq.gz/_paired.fastq.gz}
   R2unpaired=${R2//.fastq.gz/_unpaired.fastq.gz}
	
$trimmomatic_path -phred33  $R1 $R2 $R1paired $R1unpaired $R2paired $R2unpaired ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
done
