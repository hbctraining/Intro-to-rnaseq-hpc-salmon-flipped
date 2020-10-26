---
title: "Quantification of transcript abundance using Salmon"
author: "Mary Piper and Meeta Mistry"
date: "Wednesday, February 27, 2019"
---

Contributors: Mary Piper and Meeta Mistry

Approximate time: 1.25 hours

## Learning Objectives

* Explore using lightweight algorithms to quantify reads to abundance estimates
* Understand how Salmon performs quasi-mapping and transcript abundance estimation


 
## Running Salmon on multiple samples 

We just ran Salmon on a single sample, but we would like to run this on **all samples**. To do so, we will need to create a job submission script.

> *NOTE:* We are iterating over FASTQ files in the **`raw_data` directory**.

### Create a job submission script to run Salmon in serial

Since Salmon is only able to take a single file as input, one way in which we can do this is to use a 'for loop' to run Salmon on all samples in serial. What this means is that Salmon will process the dataset one sample at a time.

Let's start by opening up a script in `vim`:

```
$ vim salmon_all_samples.sbatch
```

Let's start our script with a **shebang line followed by SBATCH directives which describe the resources we are requesting from O2**. We will ask for 6 cores and take advantage of Salmon's multi-threading capabilities. 

Next we can **create a for loop to iterate over all FASTQ samples**. Inside the loop we will create a variable that stores the prefix we will use for naming output files, then we run Salmon. 

The final script is shown below:

```
#!/bin/bash

#SBATCH -p short 
#SBATCH -c 6 
#SBATCH -t 0-12:00 
#SBATCH --mem 8G 
#SBATCH --job-name salmon_in_serial 
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --reservation=HBC

cd ~/rnaseq/results/salmon

for fq in ~/rnaseq/raw_data/*.fq

do

# create a prefix
base=`basename $fq .fq`

# run salmon
salmon quant -i /n/groups/hbctraining/rna-seq_2019_02/reference_data/salmon.ensembl38.idx.09-06-2019 \
 -l A \
 -r $fq \
 -p 6 \
 -o $base.salmon \
 --seqBias \
 --useVBOpt \
 --numBootstraps 30 \
 --validateMappings

done

```

Note, that we are **adding a couple of new parameters**: 

* `-p`: specifies the number of processors or cores we would like to use for **multithreading** 
* `--numBootstraps`: specifies computation of bootstrapped abundance estimates. **Bootstraps are required for isoform level differential expression analysis for estimation of technical variance**. 
	
> _**NOTE:** `--numBootstraps` is necessary if performing **isoform-level differential expression analysis** with Sleuth, but not for gene-level differential expression analysis. Due to the statistical procedure required to assign reads to gene isoforms, in addition to the random processes underlying RNA-Seq, there will be **technical variability in the abundance estimates** output from the pseudo-alignment tool [[2](https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html), [3](https://www.nature.com/articles/nmeth.4324)] for the isoform level abundance estimates (not necessary for gene-level estimates). Therefore, **we would need technical replicates to distinguish technical variability from the biological variability** for gene isoforms._
>
> _The bootstraps estimate technical variation per gene by calculating the abundance estimates for all genes using a different sub-sample of reads during each round of bootstrapping. The variation in the abundance estimates output from each round of bootstrapping is used for the estimation of the technical variance for each gene._

Save and close the script. This is now ready to run.

```
$ sbatch salmon_all_samples.sbatch
```

---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
