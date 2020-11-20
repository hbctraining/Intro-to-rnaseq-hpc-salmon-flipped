# Workshop Schedule

> **NOTE:** The *Basic Data Skills* [Introduction to the command-line interface](https://hbctraining.github.io/Intro-to-shell-flipped/schedule/) workshop is a prerequisite.

## Pre-reading
* [Shell basics review](../lessons/shell_review.md)
* [Introduction to RNA-seq](../lessons/01_intro-to-RNAseq.md)

## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 09:45 | [Workshop introduction](../lectures/Intro_to_workshop.pdf) | Radhika |
| 09:45 - 10:25 | [Working in an HPC environment - Review](../lessons/03_working_on_HPC.md) | Radhika |
| 10:25 - 11:05 | [Project Organization and Best Practices in Data Management](../lessons/04_data_organization.md) | Meeta |
| 11:05 - 11:45 | [Quality Control of Sequence Data: Running FASTQC](../lessons/05_qc_running_fastqc_interactively.md) | Jihe |
| 11:45 - 12:00 | Overview of self-learning materials and homework submission | Jihe/Meeta |

### Self Learning #1

Before you start with the self-learning portion of the workshop, please check that **you are logged into O2** and **are working on a compute node** (i.e. your command prompt should have the word `compute` in it).

> If you are not logged into O2 or are not on a compute node, please follow the steps below as appropriate before you start with the self-learning lessons:
> 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (HBCXXcluster) (replace the "XX" in both the username and the password with the number you were assigned in class). 
> 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:00 --mem 1G /bin/bash` to get on a compute node. ***The lesson may have guidance on which arguments to modify, e.g. you may need to use more memory or more cores***.
>> ***NOTE: DO NOT specify a `reservation`.*** 
> 3. Proceed with the self learning once your command prompt has the word `compute` in it.
> 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

* [Experimental design considerations](../lessons/02_experimental_planning_considerations.md)
* [Quality Control of Sequence Data: Running FASTQC on multiple samples](../lessons/06_qc_running_fastqc_sbatch.md)
* [Quality Control of Sequence Data: Evaluating FASTQC reports](../lessons/07_qc_fastqc_assessment.md)

### Assignment #1
* Once the self-learning lessons exercises are complete, please identify the apprpopriate files as specified in the lesson (i.e. Excel file, or text/script files) and **upload** to [Dropbox](https://www.dropbox.com/request/iE82DOpC35d86Jh87zrG) the **day before the next class**.
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* Donwloadable Answer Keys: [Experimental design (one possible solution)](https://www.dropbox.com/s/524mevuyba34l5b/exp_design_table.xlsx?dl=1), [sbatch script](https://www.dropbox.com/s/9wdyhfqpic05l6p/mov10_fastqc.run?dl=1), [.out file](https://www.dropbox.com/s/l7puf8oahtbwmpk/22914006.out?dl=1), [.err file](https://www.dropbox.com/s/8a1g6o9t2kxit30/22914006.err?dl=1).

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:30 | Self-learning lessons review | All |
| 10:30 - 11:10 | [Sequence Alignment Theory](../lectures/alignment_quantification.pdf) | Meeta |
| 11:10 - 11:50 | [Quantifying expression using alignment-free methods (Salmon)](../lessons/08_quasi_alignment_salmon.md) | Radhika |
| 11:50 - 12:00 | [Review of workflow](../lectures/workflow_overview.pdf) | Radhika |

### Self Learning #2
* [Quantifying expression using alignment-free methods (Salmon on multiple samples)](../lessons/09_quasi_alignment_salmon_sbatch.md)
* QC with Alignment Data
* [Documenting Steps in the Workflow with MultiQC](../lessons/11_multiQC.md)

### Assignment #2
* Once the self-learning lessons exercises are complete, please identify the apprpopriate files as specified in the lesson (i.e. Excel file, or text/script files) and **upload** to [Dropbox](https://www.dropbox.com/request/9fWybJi6KfW8jjZXOQqK) **day before the next class**.
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
