# Workshop Schedule

> **NOTE:** The *Basic Data Skills* [Introduction to the command-line interface](https://hbctraining.github.io/Intro-to-shell-flipped/schedule/) workshop is a prerequisite.

## Pre-reading
* [Shell basics review](../lessons/shell_review.md) & Answer key
* [Introduction to RNA-seq](../lessons/01_intro-to-RNAseq.md)

## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 09:45 | Workshop introduction | Radhika |
| 09:45 - 10:25 | Working in an HPC environment | Radhika |
| 10:25 - 11:05 | Project Organization and Best Practices in Data Management | Meeta |
| 11:05 - 11:45 | Quality Control of Sequence Data: Running FASTQC | Jihe |
| 11:45 - 12:00 | Overview of self-learning materials and homework submission | Jihe/Meeta |

### Self Learning #1

Before you start with the self-learning portion of the workshop, please check that **you are logged into O2** and **are working on a compute node** (i.e. your command prompt should have the word `compute` in it).

> If you are not logged into O2 or are not on a compute node, please follow the steps below as appropriate before you start with the self-learning lessons:
> 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (HMSXXcluster) (replace the "XX" in both the username and the password with the number you were assigned in class). 
> 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:00 --mem 1G /bin/bash` to get on a compute node. ***Note, the lesson may have guidance on which arguments to modify, e.g. you may need to use more memory or more cores***.
> 3. Proceed with the self learning once your command prompt has the word `compute` in it.
> 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

* Experimental design considerations
* Quality Control of Sequence Data: Running FASTQC on multiple samples
* Quality Control of Sequence Data: Evaluating FASTQC reports

### Assignment #1
* All exercise questions from the self-learning lessons have been put together in a text file (download for local access).
  * The text file can be opened with any text editor application (i.e. Notepad++, TextWrangler) on your local computer
* Add your solutions to the exercises in the downloaded .txt file and **upload the saved text file** to Dropbox **day before the next class**.
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* Answer key

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:30 | Self-learning lessons review | All |
| 10:30 - 11:10 | Sequence Alignment Theory | Meeta |
| 11:10 - 11:50 | Quantifying expression using alignment-free methods (Salmon) | Radhika |
| 11:50 - 12:00 | Review of workflow | Radhika |

### Self Learning #2
* Quantifying expression using alignment-free methods (Salmon on multiple samples)
* QC with Alignment Data
* Documenting Steps in the Workflow with MultiQC

### Assignment #2
* All exercise questions from the self-learning lessons have been put together in a text file (download for local access).
  * The text file can be opened with any text editor application (i.e. Notepad++, TextWrangler) on your local computer
* Add your solutions to the exercises in the downloaded .txt file and **upload the saved text file** to Dropbox **day before the next class**.
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* Answer key

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:10 | Self-learning lessons review | All |
| 10:10 - 11:10 | Automating the RNA-seq workflow| Radhika |
| 11:10 - 11:45 | Troubleshooting RNA-seq Data Analysis | Meeta |
| 11:45 - 12:00 | Wrap up | Radhika |

***

## Resources



***
*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
