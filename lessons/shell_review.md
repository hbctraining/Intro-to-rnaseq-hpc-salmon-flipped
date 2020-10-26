---
title: "The Shell"
author: "Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu"
date: "October 26, 2020"
---

## Learning Objectives
- Review shell commands and concepts


## Setting up

The Introduction to RNA-seq workshop assumes that you have either a) taken our [Introduction to command-line interface workshop](https://hbctraining.github.io/Intro-to-shell-flipped/schedule/) or b) been working on the command-line and are already fluent with shell/bash. **We ask that you complete the exercises below**, to refresh some basic commands that you will be using over the course of the workshop. For each section we have relevant materials linked as a helpful reference. 

### Opening up a terminal window

*NOTE: This mandatory pre-work does not require you to login to the O2 cluster.*

On your local laptop, you will need to open up your terminal window. This will be different depending on what kind of operating system (OS) you are working on.

**With Mac OS**

Macs have a utility application called "**Terminal**" for performing tasks on the command line (shell), both locally and on remote machines. 

Please find and open the Terminal utility on your computers using the *Spotlight Search* at the top right hand corner of your screen.

**With Windows OS**

By default, there is no built-in Terminal that uses the bash shell on the Windows OS. So, we will be using a downloaded program called "**Git BASH**" which is part of the [Git for Windows](https://git-for-windows.github.io/) tool set. **Git BASH is a shell/bash emulator.** What this means is that it shows you a very similar interface to, and provides you the functionality of, the Terminal utility found on the Mac and Linux Operating systems.

Please find and open Git BASH.

> **Tip** - Windows users can use another program called [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) instead of a *bash emulator* to log in to remote machines, but it is a little more involved and has different capabilities. We encourage you to take a look at it, but we will not be covering it in this workshop.

### Downloading the example data folder

The data you will be working with can be downloaded using the link below. Clicking on the link will automatically place a file called `unix_lesson.tar.gz` to your `Downloads` folder on your computer.

- [Introduction to Shell: Dataset](https://www.dropbox.com/s/3lua2h1oo18gbug/unix_lesson.tar.gz?dl=1)

Now, in you terminal window change directories into your `Downloads` folder and check that the file is listed there:

```bash
$ cd ~/Downloads
$ ls -l unix_lesson.tar.gz 
```

To decompress the file into a folder called `unix_lesson` we use the `tar` command:

```bash
$ tar xvzf unix_lesson.tar.gz
```

- `x`: This option tells tar to extract the files.
- `v`: The “v” stands for “verbose.” This option will list all of the files one by one in the archive.
- `z`: The z option is very important and tells the tar command to uncompress the file (gzip).
- `f`: This options tells tar that you are going to give it a file name to work with (so your resulting file or folder will be named with the text before the `tar.gz` extension)

Check to see that you have the folder `unix_lesson` before proceeding.

```bash

$ ls -l unix_lesson
```

## Reviewing shell commands

### Shell basics
We are going to start this review with some basic commands pertaining to navigating around the filesystem. Helpful reference materials are listed below:

* [Introduction to Shell](https://hbctraining.github.io/Intro-to-shell-flipped//lessons/01_the_filesystem.html)
* [Wildcards and shortcuts in Shell](https://hbctraining.github.io/Intro-to-shell-flipped/lessons/02_wildcards_shortcuts.html)
* [Examining and creating files](https://hbctraining.github.io/Intro-to-shell-flipped/lessons/03_working_with_files.html)

1. Change directory into the `unix_lesson/` directory.
2. Take a quick look at the `Mov10_oe_1.subset.fq` file using `less` from `unix_lesson/`, without changing directories.
3. Move up to your home directory (parent of `unix_lesson/`).
4. With a single command change directories to the `raw_fastq/` folder.
5. With a shortest possible command change directories back to the home directory.
6. What does the `~` in the command prompt mean?
7. What is the full path to your home directory?
8. Change directories to `raw_fastq` and list all the files. (*Hint: you will need to first get into your Downloads directory*)
8. Modify the above command using the `*` wildcard to only list those files that have "oe" in their names.
10. How many and which commands have you run so far today?

### Searching and redirection
Next, we will search our files for specific patterns and redirect the results to file. Helpful reference materials are listed below:

* [Searching and redirection](https://hbctraining.github.io/Intro-to-shell-flipped/lessons/04_searching_files.html)

12. Create a new directory called `shell_review/` within the `unix_lesson/` directory.
13. Search the file `~/unix_lesson/reference_data/chr1-hg19_genes.gtf` for lines containing the string "MOV10". Save the output in the `shell_review/` directory with a new name - "Mov10_hg19.gtf".
14. Use `vim` to open the newly created file `~/unix_lesson/shell_review/Mov10_hg19.gtf` and add a comment at the top specifying how this file was created and the source of the content. Save the modified file and quit `vim`.
15. In the new file "Mov10_hg19.gtf", how many lines contain the word "exon"?

### Loops and shell scripts

* [Shell scripts and variables in Shell](https://hbctraining.github.io/Intro-to-shell-flipped/lessons/05_shell-scripts_variable.html)
* [Loops and automation](https://hbctraining.github.io/Intro-to-shell-flipped/lessons/06_loops_and_automation.html)

16. Use the `for` loop to iterate over each FASTQ file in `~/unix_lesson/raw_fastq/` and do the following:
      * Print the name of the current file
      * Generate a prefix to use for naming our output files, and store it inside a variable called `sample`.
      * Dump out the first 40 lines into a new file that will be saved in `~/unix_lesson/shell_review/`
17. Place the above `for` loop into a shell script using `vim` and run it.

**Permissions**

18. List `/n/groups/hbctraining/intro_rnaseq_hpc/` directory in long listing format
      * How many owners have files in this folder?
      * How many groups?
      * Are there any executable *files* in this folder?
      * What kind of access do you have to the `full_dataset/` directory?
      * What could user `mm573` do to take away your ability to look inside the `full_dataset/` directory?

**Environment variables**

19. Display the contents of the `$HOME` variable.
20. Use the `which` command to check where the executable file for the `pwd` command lives in the directory structure.
21. How does shell know where to find the executable file for the `pwd` command?
22. Display the contents of the variable that stores the various paths to folders containing executable command files.
23. Can you run the `bowtie2` command? What do you think you might need to do to run this command?

**LMOD system**

24. Load the `gcc/6.2.0` module.
25. Has `$PATH` changed? 
26. Load the `bowtie2/2.3.4.3` module.
27. List the modules that are loaded.

****

## Some setting up for the rest of the workshop

### Add a path to `$PATH`

We need to use one tool that is unavailable as a module on O2, but it is available in a folder on O2, so we are going to add it to our $PATH. If we just add it using the `export` command, it will only be available to us in this specific interactive session. However, if we place that export command in a script that is run everytime a new interactive session is started, it is more efficient.

* Use `vim` to open `~/.bashrc`
* Add the following line at the end of the file `export PATH=/n/app/bcbio/tools/bin:$PATH`
* Save and quit out of `vim`

### Resources on O2 and asking Slurm for them

Finally, let's review some of the information for O2 and slurm in [these slides](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lectures/HPC_intro_O2_review.pdf)

****

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
