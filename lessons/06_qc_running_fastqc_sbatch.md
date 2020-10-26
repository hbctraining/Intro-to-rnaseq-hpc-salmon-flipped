---
title: "Quality control using FASTQC"
author: "Mary Piper, Radhika Khetani"
date: Wednesday, September 20, 2017
duration: 45 minutes
---

## Learning Objectives:

* Understanding the quality values in a FASTQ file
* Using FASTQC to create a quality report
* Create and run a job submission script to automate quality assessment

## Quality Control of FASTQ files


### Performing quality assessment using job submission scripts
So far in our FASTQC analysis, we have been directly submitting commands to O2 using an interactive session (ie. `srun --pty -c 6 -p interactive -t 0-12:00 --mem 6G --reservation=HBC /bin/bash`). However, there are many more partitions available on O2 than just the interactive partition. We can submit commands or series of commands to these partitions using job submission scripts. 

**Job submission scripts** for O2 are just regular scripts, but contain the slurm **options/directives** for job submission, such as *number of cores, name of partition, runtime limit, etc*. We can submit these scripts to whichever partition we specify in the script using the `sbatch` command as follows:

```bash
# DO NOT RUN THIS
$ sbatch job_submission_script.run
```

Submission of the script using the `sbatch` command allows slurm to run your job when its your turn. Let's create a job submission script to load the FASTQC module, run FASTQC on all of our fastq files, and move the files to the appropriate directory.

Change directories to `~/rnaseq/scripts`, and create a script named `mov10_fastqc.run` using `vim`.

```bash
$ cd ~/rnaseq/scripts

$ vim mov10_fastqc.run
```

The first thing we need in our script is the **shebang line**:

```bash
#!/bin/bash
```

Following the shebang line are the O2 options. For the script to run, we need to include options for **queue/partition (-p) and runtime limit (-t)**. To specify our options, we precede the option with `#SBATCH`, which tells O2 that the line contains options for job submission to slurm. 

```bash
#SBATCH -p short 		# partition name
#SBATCH -t 0-2:00 		# hours:minutes runlimit after which job will be killed
#SBATCH -c 6 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 6G
#SBATCH --job-name rnaseq_mov10_fastqc 		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written
```
Now in the body of the script, we can include any commands we want run:

```bash
## Changing directories to where the fastq files are located
cd ~/rnaseq/raw_data

## Loading modules required for script commands
module load fastqc/0.11.3

## Running FASTQC
fastqc -t 6 *.fq

## Moving files to our results directory
mv *fastqc* ../results/fastqc/
```

Save and quit the script. Now, let's submit the job to the slurm:

```bash
$ sbatch mov10_fastqc.run
```

You can check on the status of your job with:

```bash
$ O2sacct
```

> **NOTE:** Other helpful options for checking/managing jobs are available as a [cheatsheet](https://wiki.rc.hms.harvard.edu/display/O2/O2+Command+CheatSheet) from HMS-RC.

```bash
$ ls -lh ../results/fastqc/
```
There should also be standard error (`.err`) and standard out (`.out`) files from the job listed in `~/rnaseq/scripts`. You can move these over to your `logs` directory and give them more intuitive names:

```bash
$ mv *.err ../logs/fastqc.err
$ mv *.out ../logs/fastqc.out
```

***
**Exercise**

How would you change the `mov10_fastqc.run` script if you had 9 fastq files you wanted to run in parallel.

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
