---
title: "Quality control using FASTQC - script running"
author: "Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu"
date: Friday, October 30, 2020
duration: 45 minutes
---

## Learning Objectives:

* Create and run a job submission script to automate quality assessment

## Quality Control of FASTQ files


### Performing quality assessment using job submission scripts
So far in our FASTQC analysis, we have been directly submitting commands to O2 using an interactive session (ie. `srun --pty -c 6 -p interactive -t 0-12:00 --mem 6G --reservation=HBC /bin/bash`). However, there are many [more partitions available on O2](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-Partitions(akaQueues)) than just the interactive partition. We can submit a command or series of commands to these partitions using job submission scripts. 

**Job submission scripts** for O2 are just regular shell scripts, but contain the Slurm **options/directives** for our job submission. These directives define the various resources we are requesting for our job (i.e *number of cores, name of partition, runtime limit* )

Submission of the script using the `sbatch` command allows Slurm to run your job when its your turn. Let's create a job submission script to automate what we have done in [previous lesson](05_qc_running_fastqc_interactively.md).

Our script will do the following:

1. Change directories to where the FASTQ files are located
2. Load the FastQC module
3. Run FastQC on all of our FASTQ files

Let's first change the directory to `~/rnaseq/scripts`, and create a script named `mov10_fastqc.run` using `vim`.

```bash
$ cd ~/rnaseq/scripts

$ vim mov10_fastqc.run
```

Once in the vim editor, click `i` to enter INSERT mode. The first thing we need in our script is the **shebang line**:

```bash
#!/bin/bash
```

Following the shebang line are the Slurm directives. For the script to run, we need to include options for **queue/partition (-p) and runtime limit (-t)**. To specify our options, we precede the option with `#SBATCH`. Some key resources to specify are:

|Resource|Flag|Description|
|:----:|:----:|:----:|
|partition|-p|partition name|
|time|-t|hours:minutes run limit, after which the job will be killed|
|core|-c|number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job|
|memory|--mem|memory limit per compute node for the job|

Let's specify those options as follows:

```bash
#SBATCH -p short 		# partition name
#SBATCH -t 0-2:00 		# time limit
#SBATCH -c 6 		# number of cores
#SBATCH --mem 6G   # requested memory
#SBATCH --job-name rnaseq_mov10_fastqc 		# Job name
#SBATCH -o %j.out			# File to which standard output will be written
#SBATCH -e %j.err 		# File to which standard error will be written
```

Now in the body of the script, we can include any commands we want to run. In this case, it will be the following:

```bash
## Change directories to where the fastq files are located
cd ~/rnaseq/raw_data

## Load modules required for script commands
module load fastqc/0.11.3

## Run FASTQC
fastqc -o ~/rnaseq/results/fastqc/ -t 6 *.fq
```

> **NOTE:** These are the same commands we used when running FASTQC in the interactive session. Since we are writing them in a script, the `tab` completion function will **not work**, so please make sure you don't have any typos when writing the script!

Once done with your script, click `esc` to exit the INSERT mode. Then save and quit the script by typing `:wq`. You may double check your script by typing `less mov10_fastqc.run`. If everything looks good submit the job!

```bash
$ sbatch mov10_fastqc.run
```

You should immediately see a prompt saying `Submitted batch job JobID`. Your job is assigned with that unique identifier `JobID`. You can check on the status of your job with:

```bash
$ O2sacct
```

Look for the row that corresponds to your `JobID`. The third column indicates the state of your job. Possible states include `PENDING`, `RUNNING`, `COMPLETED`. Once your job state is `RUNNING`, you should expect it to finish in less than two minutes. When the state is `COMPLETED`, that means your job is finished.

> **NOTE:** Other helpful options for checking/managing jobs are available as a [cheatsheet](https://wiki.rc.hms.harvard.edu/display/O2/O2+Command+CheatSheet) from HMS-RC.

Check out the output files in your directory:
```bash
$ ls -lh ../results/fastqc/
```
There should also be one standard error (`.err`) and one standard out (`.out`) files from the job listed in `~/rnaseq/scripts`. You can move these over to your `logs` directory and give them more intuitive names:

```bash
$ mv *.err ../logs/fastqc.err
$ mv *.out ../logs/fastqc.out
```
> **NOTE:** The `.err` and `.out` files store log information during the script running. They are helpful resources, especially when your script does not run as expected and you need to troubleshoot the script.

***
**Exercise**
1. Take a look at what's inside the `.err` and `.out` files. What do you observe? Do you remember where you see those information when using the interactive session?
2. How would you change the `mov10_fastqc.run` script if you had 9 fastq files you wanted to run in parallel? 
3. **Write out your answers in a text file and upload to the DropBox link on the schedule page.**

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
