# 01__extract_expression_data.R

# 1. Load the data
input_data_path <- "data/input/"
processed_data_path <- "data/processed/"

# CSV file of CT antigen list
cta <- data.table::fread(file = paste0(input_data_path, "ct_antigens.csv"), check.names=FALSE, stringsAsFactors=F)

# CSV file of cell surface marker list
csfm <- data.table::fread(file = paste0(input_data_path, "membrane_antigens.csv"), check.names=FALSE, stringsAsFactors=F)

# tsv.gz file of gene id and name
gene <- data.table::fread(file = paste0(input_data_path, "rna/TOIL-GTEX_TARGET_TCGA-ROW_DATA.tsv.gz"), check.names=FALSE, stringsAsFactors=F)

# TSV file of meta data
meta <- data.table::fread(file=paste0(input_data_path, "rna/TOIL-GTEX_TARGET_TCGA-COLUMN_DATA.tsv"), check.names=FALSE, stringsAsFactors=F)

# TSV file of integrated RNAseq data
rna <- data.table::fread(file=paste0(input_data_path, "rna/TOIL-GTEX_TARGET_TCGA-ASSAY_DATA-LOG2TPM.tsv"), check.names=FALSE, stringsAsFactors=F)

# 2. Data pre-processing
# 2.1. Antigen searching space
asp <- c(cta$CTantigen, csfm$Gene)
asp <- asp[asp %in% gene$gene]
asp <- unique(asp)

# 2.2. Data for further analysis
# remove data from TARGET
meta <- meta[Study!="TARGET"]
meta <- meta[Primary.site!=""]
meta <- meta[Diagnosis!=""]

# Thyroid
meta[Primary.site=="Thyroid Gland", Primary.site := "Thyroid"] 

# tissue type
meta[, type := fifelse(Diagnosis=="Normal", "normal", "cancer")]

# tissue.cancer
meta[, tissue.cancer := fifelse(Diagnosis=="Normal", Primary.site, Diagnosis)]

# remove immunoprivileged tissue
meta <- meta[Primary.site!="Testis"]

# remove column
meta[, c("Primary.site", "Diagnosis", "Gender", "Study") := NULL]

# rna data
rna <- data.table::transpose(rna, keep.names = "Sample.ID")
data.table::setnames(rna, old=colnames(rna), new=c("Sample.ID", gene$gene))

# select samples
rna <- rna[Sample.ID%in%meta$Sample.ID]

# select features
feature <- c("Sample.ID", asp)
rna <- rna[, ..feature]

rna <- merge(meta, rna, by="Sample.ID")

# 3. Result output
filename <- paste0(processed_data_path, "TCGA_GTEX_integrated_antigen_searching_space_selected.txt")
if(!exists(filename)){
    data.table::fwrite(rna, file=filename)
}
