#! /bin/bash


# This script is for PAIRED-END data which typically have the following naming convention: 
# - {sample}_R1.fastq and {sample}_R2.fastq 
# OR
# - {sample}_1.fastq and {sample}_2.fastq

# NOTE to change the extension to match the naming convention of your PE files.

for fq in ~/unix_lesson/rnaseq/raw_data/*_R1.fq
do

sbatch -p short -t 0-2:00 -c 6 --job-name rnaseq-workflow --wrap="sh ~/unix_lesson/rnaseq/scripts/PE-rnaseq_analysis_on_input_file.sh $fq"
sleep 1	# wait 1 second between each job submission
  
done
