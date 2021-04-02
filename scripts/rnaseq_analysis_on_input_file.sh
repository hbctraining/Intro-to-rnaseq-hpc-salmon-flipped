#!/bin/bash/

# This script takes a fastq file of RNA-seq data, runs FastQC, STAR, Qualimap and Salmon.
# USAGE: sh rnaseq_analysis_on_input_file.sh <name of fastq file>

# change directories to /n/scratch3/ so that all the analysis is stored there.
cd /n/scratch3/users/r/$USER/rnaseq_hbc-workshop/

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1

# grab base of filename for naming outputs
samplename=`basename $fq .subset.fq`
echo "Sample name is $samplename"

# specify the number of cores to use
cores=6

# directory with the genome and transcriptome index files + name of the gene annotation file
genome=/n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/ensembl38_STAR_index
transcriptome=/n/groups/hbctraining/rna-seq_2019_02/reference_data/salmon_index
gtf=/n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/Homo_sapiens.GRCh38.92.1.gtf

# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist
mkdir -p results/fastqc/
mkdir -p results/STAR/
mkdir -p results/qualimap/
mkdir -p results/salmon/

# set up output filenames and locations
fastqc_out=results/fastqc/
align_out=results/STAR/${samplename}
align_out_bam=results/STAR/${samplename}_Aligned.sortedByCoord.out.bam
qualimap_out=results/qualimap/${samplename}.qualimap
salmon_out=results/salmon/${samplename}.salmon
salmon_mappings=results/salmon/${samplename}_salmon.out

# set up the software environment (use version numbers)
module load fastqc/0.11.3
module load gcc/6.2.0  
module load star/2.7.0a
module load samtools/1.3.1
module load java/jdk-1.8u112
module load qualimap/2.2.1
module load salmon/1.4.0
unset DISPLAY

echo "Processing file $fq"

echo "Starting QC for $samplename"

# Run FastQC and move output to the appropriate folder
fastqc -o $fastqc_out $fq


# Run STAR
STAR --runThreadN $cores --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard

# Run Qualimap
qualimap rnaseq \
-outdir $qualimap_out \
-a proportional \
-bam $align_out_bam \
-p strand-specific-reverse \
-gtf $gtf \
--java-mem-size=8G

# Run salmon

echo "Starting Salmon run for $samplename"

salmon quant -i $transcriptome \
-p $cores \
-l A \
-r $fq \
-o $salmon_out \
--seqBias \
--useVBOpt
