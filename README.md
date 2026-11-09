# Bioinformatic_Biostatistic_Portfolio

## About me 
Master's Graduation in Bioinformatics, University of Lille. 

## Skills

Programming : R, Python, Bash, Nextflow, Snakemake, SQL, PostgrSQL  
Tools : STAR, FastQC, cutadapt, Trim Galore, CellRank2, SCENIC +, PyDESeq2, Scanpy, ScVelo, Blast  
Omics : scRNA-seq, ATAC-seq, bulk RNA-seq  
Biology : Transcriptomics, genomics, neurodevelopment  

## Projects 
[Internship : _MED13L_ multi-omics analysis 6 months](Internship_MED13L_multi_omics_analysis)

### Internship

I completed a six-month internship at the RADEME Laboratory (*Rare Diseases of Embryonic Development and Metabolism*).
I worked as a Bioinformatician within the unit under the supervision of Luc Thomès, a Bioinformatics Research Engineer.
During my internship, I worked with RNA-seq data (bulk RNA-seq and single-cell RNA-seq) as well as ATAC-seq data.

My main responsibilities involved analysing transcriptomic and epigenomic data. 
The aim was to investigate whether the MED13L mutation affected neuronal differentiation and maturation, 
potentially leading cells towards alternative differentiation pathways.

To investigate this, I analysed RNA-seq data (both single-cell and bulk RNA-seq) as well as epigenomic data (ATAC-seq).
You can find all the scripts I developed during my internship in my GitHub repository.

For the [single-cell RNA-seq](Internship_MED13L_multi_omics_analysis/Scripts/single_cell/scRNA_seq) analysis, I used CellRank2 to analyse and interpret the data. 
I developed two scripts for this purpose. 
The first one, [package_cellrank](Internship_MED13L_multi_omics_analysis/Scripts/single_cell/scRNA_seq/package_cellrank.py) contains all the functions, while the second, [run_velocity](Internship_MED13L_multi_omics_analysis/Scripts/single_cell/scRNA_seq/run_velocity.py) one allows users to provide their input data and run the functions defined in the first script.

These scripts allowed me to investigate how gene mutations can affect cell differentiation trajectories and RNA velocity.

For the [ATAC-seq](Internship_MED13L_multi_omics_analysis/Scripts/single_cell/ATAC_seq) analysis, I used SCENIC+ to investigate how chromatin accessibility and gene regulatory networks can be affected by the gene mutation.

Regarding [bulk RNA-seq](Internship_MED13L_multi_omics_analysis/Scripts/bulk_RNAseq), I chose to develop Jupyter notebooks because I was also helping a PhD student in biology develop bioinformatics skills. 
We decided to use notebooks because they are easier to understand, follow and reuse. 
This would allow her to learn how to perform the analyses and run the scripts independently after my internship.

