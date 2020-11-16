---
title: Working on HPC
author: Radhika Khetani, Meeta Mistry, Mary Piper, Jihe Liu
date: November 16, 2020
duration: 35
---

# Working in an HPC environment

## Table of contents
  * [Connect to a *login* node on O2](#connect-to-a--login--node-on-o2)
  * [Connect to a *compute* node on O2](#connect-to-a--compute--node-on-o2)
  * [More about Slurm](#more-about-slurm)
    + [Requesting resources from Slurm](#requesting-resources-from-slurm)
    + [`sbatch` job submission script](#-sbatch--job-submission-script)
  * [Using software on O2](#using-software-on-o2)
    + [LMOD system](#lmod-system)
  * [Filesystems on O2](#filesystems-on-o2)
    + [More about `/n/sctratch3`](#more-about---n-sctratch3-)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## Connect to a *login* node on O2

Type in the following command with your username to login:

```bash
ssh username@o2.hms.harvard.edu
```

You will receive a prompt for your password, and you should type in your associated password; **note that the cursor will *not move* as you type in your password**.

A warning might pop up the first time you try to connect to a remote machine, type "Yes" or "Y". 

Once logged in, you should see the O2 icon, some news, and the command prompt, e.g. `[rc_training10@login01 ~]$`.

> Note 1: `ssh` stands for secure shell. All of the information (like your password) going between your computer and the O2 login computer is encrypted when using `ssh`.

About nodes:
* A "node" on a cluster is essentially a computer in the cluster of computers. 
* There are dedicated login nodes and compute nodes.
* A login node's only function is to enable users to log in to a cluster, it is not meant to be used for any actual work/computing. There are 6 login nodes on O2.
* There are several compute nodes on O2 available for performing your analysis/work. 


## Connect to a *compute* node on O2

You can access compute nodes in 2 ways using a job scheduler or resource manager like Slurm.
1. Directly using an interactive session (Slurm's `srun` command): 
    * The `srun` command with a few mandatory parameters will create an "interactive session" on O2. 
    * This is essentially a way for us to do work on the compute node directly from the terminal. 
    * If the connectivity to the cluster is lost in the middle of a command being run that work will be lost in an interactive session.
1. By starting a "batch" job (Slurm's `sbatch` command): 
    * The `sbatch` command with a few mandatory parameters + a specialized shell script will result in the script being run on a compute node. 
    * This "job" will not be accessible directly from the Terminal and will run in the background. 
    * Users do not need to remain connected to the cluster when such a "batch job" is running.

For now let's start an interactive session on O2 using `srun`. 

```bash
$ srun --pty -p interactive -t 0-3:00 --mem 1G --reservation=HBC1 /bin/bash
```

In the above command the parameters we are using are requesting specific resources:
* `--pty` - Start an interactive session
* `-p interactive` - on the "partition" called "interactive" (a partition is a group of computers dedicated to certain types of jobs, interactive, long, short, high-memory, etc.)
* `-t 0-8:00` - time needed for this work: 0 days, 8 hours, 0 minutes.
* `--mem 1G` - memory needed - 1 gigabyte
* `--reservation=HBC1` - *this is only for **in class** portions of this workshop, make sure you don't use it for self-learning or when you have your own accounts.*
* `/bin/bash` - You want to interact with the compute node using the *bash* shell

> These parameters are used for `sbatch` as well, but they are listed differently within the script used to submit a batch job. We will be reviewing this later in this lesson.

**Make sure that your command prompt now contains the word "compute", *e.g. `[rc_training10@compute-a-16-163 ~]$`*.** 
You are now working on a compute node directly in an "interactive" session!

Let's check how many jobs we have running currently, and what resources they are using.

```bash
$ O2squeue
```

## More about Slurm

* Slurm = Simple Linux Utility for Resource Management
* Fairly allocates access to resources (computer nodes) to users for some duration of time so they can perform work
* Provides a framework for starting, executing, and monitoring batch jobs
* Manages a queue of pending jobs; ensures that no single user or core monopolizes the cluster

### Requesting resources from Slurm

Below is table with some of the arguments you can specify when requesting resources from Slurm for both `srun` and `sbatch`:

| Argument | Description / Input | Examples | Links |
|:-----------:|:----------:|:--------:|:----------:|
| -p | name of compute partition | short, medium, interactive | [O2 Wiki - Guidelines for choosing a partition](https://wiki.rc.hms.harvard.edu/display/O2/How+to+choose+a+partition+in+O2) | 
| -t | how much time to allocate to job | 0-03:00, 5:00:00 | [O2 Wiki - Time limits for each partition](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-Timelimits) |
| -c | max cores | 4, 8, 20 | [O2 Wiki - How many cores?](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-Howmanycores?) |
| --mem | max memory | 8G, 8000 | [O2 Wiki - Memory requirements](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-Memoryrequirements) |
| -o | name of file to create with standard output | %j.out | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| -e | name of file to create with standard error | %j.err | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| -J | name of the job | Fastqc_run, rnaseq_workflow_mov10 | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| --mail-type | send an email when job starts, ends or errors out  | END, ALL | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| --mail-user | send email to this address | xyz10@harvard.edu | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |

### `sbatch` job submission script

An `sbatch` job submission script is essentially a normal shell script with the Slurm resources specified at the top (Slurm directives). Below is an example of an sbatch shell script that is requesting the following: 
* the "short" partition for 2 hours 
* on 4 cores (30 minutes for each core)
* using 400MBs (100MB for each core)

***DO NOT RUN***
```
#! /bin/sh

#SBATCH -p short
#SBATCH –t 0-02:00
#SBATCH –c 4
#SBATCH --mem=400M
#SBATCH –o %j.out
#SBATCH –e %j.err
#SBATCH -J fastqc_run
#SBATCH --mail-type=ALL
#SBATCH –-mail-user=xyz10@med.harvard.edu

## Load the fastqc module
module load fastqc/0.11.5

# run fastqc (multithreaded)
fastqc -t 4 file1_1.fq file1_2.fq file2_1.fq file2_2.fq
```

## Using software on O2

### LMOD system

In the above example we have used the FastQC tool, and prior to using it we have used the command `module load fastqc/0.11.5`. The `module load` command is part of the LMOD system available on O2. It enables users to access software installed on O2 easily, and manages every software's dependency. LMOD system adds directory paths of software and their dependencies (if any) into the `$PATH` variable.

So, instead of using `/n/app/fastqc/0.11.5/bin/fastqc` as our command, we can load the module and use `fastqc` as the command. 

Some key LMOD commands are listed below:

| LMOD command | description |
|:---------:|:---------:|
| `module spider` | List all possible modules on the cluster |
| `module spider modulename` | List all possible versions of that module |
| `module avail` | List available modules available on the cluster |
| `module avail string` | List available modules containing that string |
| `module load modulename/version` | Add the full path to the tool to `$PATH` |
| `module list` | List loaded modules | 
| `module unload modulename/version` | Unload a specific module |
| `module purge` | Unload all loaded modules |

> Note: On O2, a lot of tools used for analysis of sequencing data need to have the `gcc` compiler module loaded (`module load gcc/6.2.0`) prior to loading the tool of interest.

*** 
**Exercise**

1. What are the contents of the `$PATH` environment variable?
1. Try running the `multiqc` command. What do you get as output?
1. Check if the `multiqc` tool is available as a module. How many versions of `multiqc` are available?
1. Load `multiqc/1.9`. Did you have to load any additional modules?
1. List all the modules loaded in your environment
1. Try running the `multiqc` command again. What do you get as output?
1. Check the contents of the `$PATH` environment variable again. Any changes from before?

***

## Filesystems on O2

<p align="center">
<img src="../img/O2_primary-storage.png" width="600">
</p>

* Storage on HPC systems is organized differently than on your personal machine.
* Each node on the cluster does not have storage; instead, it is on disks bundled together externally.
* Storage filesystems can be quite complex, with large spaces dedicated to a pre-defined purpose.
* Filesystems are accessed over the internal network by all the nodes on the cluster.
* There are 3 major groups on the cluster, each with their features and constraints:
   1. `/n/data1`, `/n/data2`, `/n/groups` - Large datasets are stored in these parent directories (see features/constraints in the image above).
   2. `/home` - the home directories of all users are under this parent directory (see features/constraints in the image above).
   3. `/n/sctratch3` - scratch space for temporary storage.

### More about `/n/sctratch3`

* It is for data only needed temporarily during analyses.
* Each user can use up to 10 TB and 1 million files/directories.
* Files not accessed for 30 days are automatically deleted.
* **No backups!**
* You can create your own folder using this command `/n/cluster/bin/scratch3_create.sh`

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
