UCSC Xena (TOIL) re-processed RNAseq data from GTEX (normal tissues), TARGET and TCGA (tumour samples only)

## File descriptions

### TOIL-GTEX_TARGET_TCGA-COLUMN_DATA.tsv.gz

Column (sample) data. 18392 rows (excluding header), 5 columns:

	**Sample ID** - unique sample identifier, same as the column names in the assay file
	**Primary site** - anatomical origin of the sample
	**Diagnosis** - either 'Normal' or tumour type
	**Gender** - donor's gender ('female' or 'male' only)
	**Study** - data source ('GTEX', 'TARGET' or 'TCGA')

### TOIL-GTEX_TARGET_TCGA-ROW_DATA.tsv.gz

Row (gene) data. 60498 rows (excluding header), 6 columns:

	**id** - ENSEMBL gene identifier, with point version. Order is the same as in the assay file
	**gene** - gene symbol (HUGO?)
	**chrom** - chromosome (no alternative contigs)
	**chromStart** - gene's upstream (left-most) position on the chromosome
	**chromEnd** - gene's downstream (right-most) position on the chromosome
	**strand** - gene's orientation ('+' or '-')

### TOIL-GTEX_TARGET_TCGA-ASSAY_DATA-LOG2TPM.tsv.gz

Assay data in log2(TPM + 0.001) units. Gene order is the same as in the `id` column of the gene table. Column names are sample IDs (same as in `Sample ID` column of sample data).
60498 rows, 18392 columns (excluding header).

### TOIL-GTEX_TARGET_TCGA-SE.Rds

This is an R object of SummarizedExperiment class, contains assay data, as well as column (as DFrame) and row data (as GRanges object).
