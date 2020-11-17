---
title: "Quality control using FASTQC"
author: "Mary Piper, Radhika Khetani, Meeta Mistry, Jihe Liu"
date: Friday, October 30, 2020
duration: 45 minutes
---

## Learning Objectives:

* Understand the quality values in a FASTQ file
* Create a quality report using FASTQC
 
## Quality Control of FASTQ files


The first step in the RNA-Seq workflow is to take the FASTQ files received from the sequencing facility and assess the quality of the sequence reads. 

<p align="center">
<img src="../img/rnaseq_salmon_workflow.png" width="400">
</p>

### Unmapped read data (FASTQ)

The [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file format is the defacto file format for sequence reads generated from next-generation sequencing technologies. This file format evolved from FASTA in that it contains sequence data, but also contains quality information. Similar to FASTA, the FASTQ file begins with a header line. The difference is that the FASTQ header is denoted by a `@` character. For a single record (sequence read), there are four lines, each of which are described below:

|Line|Description|
|----|-----------|
|1|Always begins with '@', followed by information about the read|
|2|The actual DNA sequence|
|3|Always begins with a '+', and sometimes the same info as in line 1|
|4|Has a string of characters representing the quality scores; must have same number of characters as line 2|

Let's use the following read as an example:

```
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
```

The line 4 has characters encoding the quality of each nucleotide in the read. The legend below provides the mapping of quality scores (Phred-33) to the quality encoding characters. *Different quality encoding scales exist (differing by offset in the ASCII table), but note the most commonly used one is fastqsanger, which is the scale output by Illumina since mid-2011.* 
 ```
 Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                   |         |         |         |         |
    Quality score: 0........10........20........30........40                                
```
 
Using the quality encoding character legend, the first nucelotide in the read (C) is called with a quality score of 31 (corresponding to encoding character `@`), and our Ns are called with a score of 2 (corresponding to encoding character `#`). **As you can tell by now, this is a bad read.** 

Each quality score represents the probability that the corresponding nucleotide call is incorrect. This quality score is logarithmically based and is calculated as:

	Q = -10 x log10(P), where P is the probability that a base call is erroneous

These probabaility values are the results from the base calling algorithm and dependent on how much signal was captured for the base incorporation. The score values can be interpreted as follows:

|Phred Quality Score |Probability of incorrect base call |Base call accuracy|
|:-------------------:|:---------------------------------:|:-----------------:|
|10	|1 in 10 |	90%|
|20	|1 in 100|	99%|
|30	|1 in 1000|	99.9%|
|40	|1 in 10,000|	99.99%|

Therefore, for the first nucleotide in the read (C), there is less than a 1 in 1000 chance that the base was called incorrectly. Whereas, for the the end of the read there is greater than 50% probabaility that the base is called incorrectly.

## Assessing quality with FastQC

Now that we understand what information is stored in a FASTQ file, the next step is to examine quality metrics for our data.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides a simple way to do some quality checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses, which you can use to obtain an impression of whether your data has any problems that you should be aware of before moving on to the next analysis.

FastQC does the following:
* accepts FASTQ files (or BAM files) as input
* generates summary graphs and tables to help assess your data
* generates an easy-to-view HTML-based report with the graphs and tables

***

> NOTE: Before we run FastQC, **you should be on a compute node** in an interactive session. Please run the following `srun` command if you are not on a compute node.
> 
> ```bash
> $ srun --pty -p interactive -t 0-3:00 --mem 1G --reservation=HBC1 /bin/bash
> ```
>
> ***An interactive session is very useful to test tools and workflows.***

### Run FastQC  

Change directories to `raw_data`.

```bash
$ cd ~/rnaseq/raw_data
```  

Before we start using software, we have to load the module for each tool. On O2, this is done using an **LMOD** system. 

If we check which modules we currently have loaded, we should not see FastQC.

```bash
$ module list
```

This is because the FastQC program is not in our $PATH (i.e. it's not in a directory that shell will automatically check to run commands/programs).

```bash
$ echo $PATH
```

To run the FastQC program, we first need to load the appropriate module, so it puts the program into our path. To find the FastQC module to load we need to search the versions available:

```bash
$ module spider fastqc
```

Once we know which version we want to use (0.11.3), we can load the FastQC module:

```bash
$ module load fastqc/0.11.3
```

Once a module for a tool is loaded, you have essentially made it directly available to you like any other basic shell command.

```bash
$ module list

$ echo $PATH
```

Now, let's create a directory to store the output of FastQC:

```bash
$ mkdir ~/rnaseq/results/fastqc
```

We will need to specify this directory in the command to run FastQC. How do we know which argument to use?

```bash
$ fastqc --help
```

> **NOTE:** From the help manual, we know that `-o` (or `--outdir`) will create all output files in the specified output directory. Note that another argument, `-t`, specifies the number of files which can be processed simultaneously. We will use `-t` argument later. You may explore other arguments as well based on your needs.

FastQC will accept multiple file names as input, so we can use the `*.fq` wildcard.

```bash
$ fastqc -o ~/rnaseq/results/fastqc/ *.fq
```

*Did you notice how each file was processed serially? How do we speed this up?*

FastQC has the capability of splitting up a single process to run on multiple cores! To do this, we will need to specify an additional argument `-t` indicating number of cores. We will also need to exit the current interactive session, since we started this interactive session with only 1 core. We cannot have a tool to use more cores than requested on a compute node. 

Exit the interactive session and start a new one with 6 cores:

```bash
$ exit  #exit the current interactive session (you will be back on a login node)

$ srun --pty -c 6 -p interactive -t 0-3:00 --mem 2G --reservation=HBC1 /bin/bash  #start a new one with 6 cores (-c 6) and 2GB RAM (--mem 2G)
```

Once you are on the compute node, check what job(s) you have running and what resources you are using.

```bash
$ O2squeue
```

Now that we are in a new interactive session with the appropriate resources, we will need to load the module again for this new session.

```bash
$ module load fastqc/0.11.3  #reload the module for the new (6-core) interactive session
```

We will also move into the `raw_data` directory (remember we are on a new compute node now):

```bash
$ cd ~/rnaseq/raw_data
```

Run FastQC and use the multi-threading functionality of FastQC to run 6 jobs at once (with an additional argument `-t`).

```bash
$ fastqc -o ~/rnaseq/results/fastqc/ -t 6 *.fq  #note the extra parameter we specified for 6 threads
```

*Do you notice a difference? Is there anything in the ouput that suggests this is no longer running serially?*

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
