# Working in an HPC environment


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

## More about Slurm

* Slurm = Simple Linux Utility for Resource Management
* Fairly allocates access to resources (computer nodes) to users for some duration of time so they can perform work
* Provides a framework for starting, executing, and monitoring batch jobs
* Manages a queue of pending jobs; ensures that no single user or core monopolizes the cluster

### Requesting resources from Slurm

Below is table with some of the arguments you can specify when requesting resources from Slurm for both `srun` and `sbatch`:

| Argument | Description / Input | Examples | Links |
|:-----------:|:----------:|:--------:|:----------:|
| -p | name of compute partition | short | [O2 Wiki - Guidelines for choosing a partition](https://wiki.rc.hms.harvard.edu/display/O2/How+to+choose+a+partition+in+O2) | 
| -t | how much time to allocate to job | 0-03:00 | [O2 Wiki - Time limits for each partition](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-Timelimits) |
| -c | max cores | 4 | [O2 Wiki - How many cores?](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-Howmanycores?) |
| --mem | max memory | 8G | [O2 Wiki - Memory requirements](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-Memoryrequirements) |
| -o | name of file to create with standard output | %j.out | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| -e | name of file to create with standard error | %j.err | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| -J | name of the job | bowtie2_run1 | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| --mail-type | send an email when job starts, ends or errors out  | END | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |
| --mail-user | send email to this address | xyz10@harvard.edu | [O2 Wiki - `sbatch` options quick reference](https://wiki.rc.hms.harvard.edu/display/O2/Using+Slurm+Basic#UsingSlurmBasic-sbatchoptionsquickreference) |

### `sbatch` job submission script



## Using software on O2

### LMOD system

* Load the `gcc/6.2.0` module.
* Has `$PATH` changed? 
* Load the `bowtie2/2.3.4.3` module.
* List the modules that are loaded.

## Filesystems on O2

### Scratch space

We need to use one tool that is unavailable as a module on O2, but it is available in a folder on O2, so we are going to add it to our $PATH. If we just add it using the `export` command, it will only be available to us in this specific interactive session. However, if we place that export command in a script that is run everytime a new interactive session is started, it is more efficient.

* Use `vim` to open `~/.bashrc`
* Add the following line at the end of the file `export PATH=/n/app/bcbio/tools/bin:$PATH`
* Save and quit out of `vim`

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
