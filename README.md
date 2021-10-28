## Orthogonal CRISPR-Cas tools for genome editing, inhibition, and CRISPR recording in zebrafish embryos 
### Paige R. Takasugi, Shengzhou Wang, Kimberly T. Truong, Evan P. Drage, Sahar N. Kanishka, Marissa A. Higbee, Nathan Bamidele, Ogooluwa Ojelabi, Erik J. Sontheimer, James A. Gagnon 

This repository contains all the necessary code and data files to reproduce the following figure, Figure S4, and S5 of the manuscript. 

Below is Figure 4, which exemplifies the potential of our orthogonal CRISPR system for in vivo recording. Parts C-F can be created using this repository. 

**The basic outline is as follows:**
  1. Identify high-quality CRISPR-induced mutations in barcodes by running the GESTALT pipeline (McKenna et al. 2016).
  2. Process outputs of edited barcodes to generate main and supplemental figures. 

**If you want to generate the figures without spending a lot of time and compute for step 1, skip to section [For generating figures using processed data](https://github.com/Gagnon-lab/takasugi-genetics/#For-generating-figures-using-processed-data)** 

## Data pre-processing using GESTALT pipeline 
To call CRISPR-induced mutations, we have adapted the pipeline from https://github.com/aaronmck/Cas9FateMapping. Briefly, this takes raw sequencing data, merges paired reads, aligns resulting consensus sequences to the reference, and generates output files for further processing and downstream analysis. 

To generate processed outputs from sequencing data, download the raw FASTQ files from GEO accession number [GSE186338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186338).

We recommend running this pipeline on a high-performance computing cluster as it will require quite a bit of memory (and time) to execute. 

Although most of the instructions can be found [here](https://github.com/mckennalab/SingleCellLineage), we briefly outline the instructions to get started on a HPC. 

# Dependencies
* Docker-like container 
* Javascript/HTML/D3

First, download the Docker container. The Center for High Performance Computing (CHPC) at the University of Utah offers Singularity as a container environment and one can download the container as a Singularity image like so: 

```
singularity build GESTALT.sif docker://aaronmck/genomics:sc_GESTALT 
```

The rest of the files needed to run the pipeline can be mostly found in `pipeline/` and are as follows:
1. `barcode.referenceseq.fa`: fasta file containing the barcode sequence including
2. `barcode.referenceseq.fa.primers`: primer sequences expected to be on both ends of the amplicon for sequencing, one on each line. 
3. `barcode.referenceseq.fa.cutSites`: tab-delimited file with three columns; the first indicates the sequence of the CRISPR target, the second is the start position (5') of that target, and the third is the position of the predicted cutsite depending on the Cas system. 
4. `gestalt_pipeline_no_trees.scala`: pipeline script that is dependent on the GATK Queue processing engine 
5. `embryo-samples.tearsheet.txt`: tearsheet listing all the samples to create compute jobs for the Queue manager 
6. `run_gestalt_pipeline.sh`: shell script to execute the pipeline with input samples provided in the tearsheet
7. Depending on whatever is used as a job management system on your cluster, one needs to provide a SLURM-like script that at a minimum does the following: 
```
# makes the output folder to dump outputs
mkdir -p OUTPUT

# executes the pipeline after mounting the data/pipeline files into the container 
singularity -exec --bind $DATA_DIR:/my_data GESTALT.sif /my_data/run_gestalt_pipeline.sh

```

One should copy all the above files and FASTQs into `$DATA_DIR`. 

## For generating figures using processed data

