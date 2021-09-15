---
title: "Sequence Alignment Theory"
author: "Radhika Khetani, Meeta Mistry, Will Gammerdinger"
date: Tuesday, October 12, 2021
duration: XX minutes
---

## Learning Objectives:

* Describe differences between genome-mapping versus transcriptome-mapping techniques when aligning RNA-seq data
* Differentiate between the various reference files needed when running RNA-seq alignment

## What is Sequence Alignment?

Sequence alignment is a process through which DNA, RNA or protein sequences are compared to a reference sequence(s) in order to find a corresponding match between the queried sequence and the reference sequence. You may have already come across the popular **BLAST** (<ins>B</ins>asic <ins>L</ins>ocal <ins>A</ins>lignment <ins>S</ins>earch <ins>T</ins>ool) algorithm already. Sequence alignment is one of the most fundamental tools in bioinformatics and it has wide-ranging applications including:

  * Deriving the likely species of origin for an unknown sample
  * Comparing genomic rearrangements between closely related species
  * Estimating gene expression
  * Many more applications
 
 In addition to [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), other popular alignment tools include [HISAT2](http://daehwankimlab.github.io/hisat2/), [bwa](http://bio-bwa.sourceforge.net), [STAR](https://github.com/alexdobin/STAR) and others. Each of these alignement tools handles sequence alignment differently. Some of these differences are important and we will discuss below, but many of these differences are beyond the scope of this lesson. Recently, some software packages used for analyzing RNA-Seq data ([Salmon](https://combine-lab.github.io/salmon/) and [Kallisto](https://pachterlab.github.io/kallisto/about)) have been introduced and they handle alignment a bit differently than the previously mentioned alignment tools. In this lesson, we will discuss historical and current alignment techniques in order to prepare you for carrying out your own sequence alignment.
 
 



