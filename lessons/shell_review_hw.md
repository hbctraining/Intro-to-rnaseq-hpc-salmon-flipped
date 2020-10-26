---
title: "The Shell"
author: "Mary Piper, Radhika Khetani, Meeta Mistry"
date: "February 2019"
---

## Learning Objectives
- Review shell commands and concepts
- Review HPC terminology and commands

## Setting up

**With Macs**

Macs have a utility application called "**Terminal**" for performing tasks on the command line (shell). We can open this utility to access the shell.

**With Windows**

By default, there is no terminal for the bash shell available in the Windows OS, so you have to use a downloaded program, "**Git BASH**". Git BASH is part of the [Git for Windows](https://git-for-windows.github.io/) download, and is a shell (bash) emulator.

#### Command prompt

Once you have opened the shell, you should see the command prompt ending with `$`. It will have some characters before the `$`, something like `[MacBook-Pro-5:~]`, this is telling you what the name of the computer you are working on is. 

```bash
[MacBook-Pro-5:~]$ 
```

### Downloading data

In order to go through the exercises and review some commands and concepts, we will be working with some RNA-Seq data. We need to **download the data to our current folder** using the link below. To know what folder we are currently inside, we can use the 'print working directory' command:

```bash
$ pwd
```

On a **Mac** your current folder should be something starting with `/Users/`, like `/Users/marypiper/`.

On a **Windows** machine, our current folder should be something starting with `C:/Users/marypiper`. To find this in your File explorer try clicking on PC and navigating to that path.

Once you have identified which folder you are in, this is where we will be downloading your data. Right click on the link below, and be sure to **Save link as**. You will be prompted to choose a folder within a Finder/File Explorer window. Navigate to the directory that was listed when running `pwd`.

**Download RNA-Seq data to your working directory:** right-click [here](https://github.com/hbctraining/Training-modules/blob/master/Intro_shell/data/unix_lesson.zip?raw=true) and choose **Save link as**.

If you have downloaded the file to the correct location, type the 'list' command:

```bash
$ ls
```

You should see `unix_lesson.zip` as part of the output to the screen.

Finally, to decompress the folder, we can use the `unzip` command:

```bash
$ unzip unix_lesson.zip 
```

You should see output stating the contents of the folder are being decompressed or inflated; this is good. Now when you run the `ls` command again you should see a folder called `unix_lesson`.

```bash
$ ls
```

## Reviewing shell commands

Below are some questions that you can work on to review your knowledge of bash. Remember, there are likely multiple ways to do the same thing!

****

**Shell basics**

1. Change directory into the `unix_lesson/` directory.
2. Take a quick look at the `Mov10_oe_1.subset.fq` file using `less` from `unix_lesson/`, without changing directories.
3. Move up to your home directory (parent of `unix_lesson/`).
4. With a single command change directories to the `raw_fastq/` folder.
5. With a shortest possible command change directories back to the home directory.
6. What does the `~` in the command prompt mean?
7. What is the full path to your home directory?
8. Use the `*` wildcard to only list those files in `raw_fastq/` that have "oe" in their names.
9. List the last 10 commands you have run so far in this lesson.

**Searching and redirection**

10. Move back in to the `unix_lesson/` directory. Create a new directory called `shell_review/` within the `unix_lesson/` directory.
11. Grab the lines in `reference_data/chr1-hg19_genes.gtf` with the string "MOV10" in them and save the output in the `shell_review/` directory with a new name - "Mov10_hg19.gtf".
12. Use `vim` to open the newly created file `unix_lesson/shell_review/Mov10_hg19.gtf` and add a comment at the top specifying how this file was created and the source of the content. Save the modified file and quit `vim`.
13. How many lines in the new file have the word "exon" in them?

**Loops and shell scripts**

14. Use the `for` loop to iterate over each FASTQ file in `raw_fastq/` and do the following:
      * Generate a prefix to use for naming our output files
      * Print the name of the current file
      * Dump out the first 40 lines into a new file that will be saved in `~/unix_lesson/shell_review/`
15. Place the above `for` loop into a shell script using `vim` and run it.


**Permissions** 

16. You are working on the cluster and are accessing files in the `/n/groups/hbctraining/` directory in long listing format:

```bash
$ ls -l /n/groups/hbctraining/

drwxrwsr-x  4 mm573 hbctraining 831 Feb 29  2016 bcbio-rnaseq
drwxrwsr-x 12 mm573 hbctraining 318 May 24 11:13 chip-seq
-rw-r--r--  1 root  hbctraining   0 Apr  5  2015 copy_me.txt
drwxrwsr-x  3 rsk27 hbctraining 201 Apr  5  2015 exercises
drwxrwsr-x  6 rsk27 hbctraining 293 Oct 27 09:40 for_chipseq
drwxrwsr-x  4 mp298 hbctraining  51 Dec  6  2016 mep-data
drwxrwsr-x  4 rsk27 hbctraining  53 Jun  2 15:57 ngs_course
drwxrwsr-x  4 rsk27 hbctraining  53 Nov  2  2016 ngs-course_backup_Nov1-2016
drwxrwsr-x  6 mm573 hbctraining 107 Mar 24  2016 ngs-data-analysis2016
```

Answer the following questions based on information provided above:

* How many owners have files in this folder?
* How many groups?
* Are there any executable *files* in this folder?
* What could user `mm573` do to take away your ability to look inside the `chip-seq/` directory?

**Environment variables** 

17. What is the importance of the `$HOME` variable?
18. How does shell know where to find the executable file for the `pwd` command?
19. Display the contents of the variable that stores the various paths to folders containing executable command files.

****

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
