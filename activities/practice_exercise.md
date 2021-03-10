_**To perform this exercise you will need an O2 account. You can request an account by following the instructions on [O2's account request page](https://harvardmed.service-now.com/stat?id=service_catalog_cards&sys_id=5165e1dbdb209050b642f27139961979&sysparm_category=991a7f2edb890c10b642f2713996196a).**_

## Running the RNA-seq workflow

We have downloaded the raw FASTQ files from the SRA for the sequencing data used in the paper: [Silencing SMOC2 ameliorates kidney fibrosis by inhibiting fibroblast to myofibroblast transformation](https://pubmed.ncbi.nlm.nih.gov/28422762/). The paper explores kidney fibrosis in wildtype and SMOC2-overexpressing mice. 

>_**NOTE:** If you are interested in downloading other datasets from the SRA, we have [materials](https://hbctraining.github.io/Accessing_public_genomic_data/lessons/downloading_from_SRA.html) available detailing how to do this._

### Set-up
1. Copy the compressed experimental data folder from `/n/groups/hbctraining/kidney_fibrosis_rnaseq.tar.gz` to your own `/n/scratch3/users/ecommonsID` directory.
2. Extract the directory using the command `tar -xzvf kidney_fibrosis_rnaseq.tar.gz`. This command may take a while to run.
3. Look inside the directory, you should find the following:

    - a `raw_fastq` folder containing the raw fastq files
    - a `meta` folder with a metadata file containing information about each of the samples
4. Create a `reference_data` folder and download the transcriptome FASTA file for mouse to the folder. 

    - For Ensembl references, go to [http://useast.ensembl.org/info/data/ftp/index.html](http://useast.ensembl.org/info/data/ftp/index.html)
    - Find the mouse species row and click on the *FASTA* link in the **cDNA (FASTA)** column. 
    - Right-click on the link for the `*cdna.all.fa.gz` file to copy it.
    - Navigate to the `reference_data` folder and run the command `wget <paste contents of link>`. This should download the transcriptome FASTA file to the directory.
    - Extract the `*cdna.all.fa.gz` file by running the code: `gzip -d *cdna.all.fa.gz`.
5. Set-up additional expected folders (e.g. results, etc.) for your project (i.e. create subdirectories and additional directories where you feel is necessary). 

### Analysis
Using the workflow and submission scripts we generated in class, parallelize the RNA-Seq analysis of all files in this dataset. For each FASTQ file you will need to perform the following:

  - Run FastQC
  - Generate abundance estimates with Salmon
  - Evaluate the MultiQC report
 
  **HINT: You will need to create a mouse index for Salmon.** 
