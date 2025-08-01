## Method:
- Antigen search space: a list of 1610 genes with a cancer-testis (CT) expression pattern was assembled from different resources (manually curated literature, text mining and high-throughput screens) as the candidate antigen for SPEAR-T cells. The potential target antigens for TRuC-T cells are defined as genes with known or predicted cell surface expression from the COMPARTMENTS database 
- Gene expression data: RNA-seq data across 9795 samples taken from TCGA and 7425 samples from GTEx were extracted in the TPM format from the UCSC database, which had been uniformly processed by the Toil pipeline
- Intelligent subsampling partitioning was used to increase the speed of antigen identification. Samples were partitioned using geometric sketching to get an equal representation of all tissue types and tumour samples in both partitions
- Clustering-based method for antigen pair identification Davies-Bouldin index, which measures the ratio of within-cluster spread to between-cluster distance, was employed to quantify the separation between samples of one type of tumour versus all normal tissues in training dataset.

<img width="326" height="133" alt="image" src="https://github.com/user-attachments/assets/3ace9c9a-7d89-487f-9ee8-97f119862ef0" />

#

<img width="175" height="172" alt="image" src="https://github.com/user-attachments/assets/cdcca350-5eb6-47a1-9d00-035c39ab0384" />
