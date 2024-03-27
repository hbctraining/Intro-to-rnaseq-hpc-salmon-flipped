#!/bin/bash

#SBATCH -p short 
#SBATCH -c 6 
#SBATCH -t 0-12:00 
#SBATCH --mem 8G 
#SBATCH --job-name salmon_in_serial 
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=xyz10@harvard.edu

# Load Salmon module
module load salmon/1.8.0

# Change directory to where the Salmon results will be output
cd ~/rnaseq/results/salmon

# Main script for running salmon with for loop

for fq in ~/rnaseq/raw_data/*.fq

do

# create a prefix for the output file
samplename=`basename $fq .fq`

# run salmon
salmon quant -i /n/groups/hbctraining/rna-seq_2019_02/reference_data/salmon_index \
 -l A \
 -r $fq \
 -o ${samplename}_salmon \
 --seqBias \
 --useVBOpt \
 --validateMappings \
 -p 6 \

done

