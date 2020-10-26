#!/bin/bash

#SBATCH -p priority             # partition name
#SBATCH -t 0-12:00              # hours:minutes runlimit after which job will be killed
#SBATCH -c 6            # number of cores requested - what you plan to use to run your job
#SBATCH --mem 8G
#SBATCH --job-name salmon_mapping_PE           # Job name
#SBATCH -o salmon-mapping.out                       # File to which standard out will be written
#SBATCH -e salmon_mapping.err               # File to which standard err will be written


# Change directories
cd cd ~/rnaseq/raw_data

# Get all sample names from a file that contains the prefix
files=`cut -f 1 samples.csv`

for sample in $files

  do

  salmon quant -i /n/groups/hbctraining/rna-seq_2019_02/reference_data/salmon.ensembl38.idx \
     -l A \
     -r ${sample}_R1.fastq ${sample}_R2.fastq \
     -o ../results/salmon/${sample} \
     -p 6 \
     --seqBias \
     --useVBOpt \
     --numBootstraps 30

done
