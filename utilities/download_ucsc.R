## Download the data

### TCGA
dir.create("TCGA")
if (!file.exists("TCGA/tcga_target_no_normal_rsem_gene_tpm.gz")) {
      download.file(
            "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_target_no_normal_rsem_gene_tpm.gz",
            destfile = "TCGA/tcga_target_no_normal_rsem_gene_tpm.gz")
}

if (!file.exists("TCGA/TCGA_TARGET_phenotype")) {
      download.file(
            "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_TARGET_phenotype",
            destfile = "TCGA/TCGA_TARGET_phenotype")
}


if (!file.exists("TCGA/gencode.v23.annotation.gene.probemap")) {
      download.file(
            "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap",
            destfile = "TCGA/gencode.v23.annotation.gene.probemap")
}

### GTEX
dir.create("GTEX")
if (!file.exists("GTEX/gtex_RSEM_gene_tpm.gz")) {
      download.file(
            "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_gene_tpm.gz",
            destfile = "GTEX/gtex_RSEM_gene_tpm.gz")
}

if (!file.exists("GTEX/GTEX_phenotype.gz")) {
      download.file(
            "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/GTEX_phenotype.gz",
            destfile = "GTEX/GTEX_phenotype.gz")
}


## Process and combine colData - clinical/sample info
colData <- data.frame(
   `Sample ID` = character(),
   `Primary site` = character(),
   `Diagnosis` = character(),
   `Gender` = character(),
   `Study` = character(),
   check.names = F, stringsAsFactors = F
)

### TCGA/TARGET
TT <- read.table("TCGA/TCGA_TARGET_phenotype", sep="\t", header = T, check.names = F, stringsAsFactors = F)
TT <- TT[, c("sampleID", "_primary_site", "_primary_disease", "gender", "_study")]
names(TT) <- names(colData)
TT$Gender <- tolower(TT$Gender)
TT$Study <- ifelse(grepl("^TCGA", TT$`Sample ID`), "TCGA", ifelse(
   grepl("^TARGET", TT$`Sample ID`), "TARGET", TT$Study
))

### GTEX
GT <- read.table("GTEX/GTEX_phenotype.gz", sep="\t", header = T, check.names = F, stringsAsFactors = F)
GT <- GT[!grepl("Cells", GT$`body_site_detail (SMTSD)`), ]  # Remove cell line data
GT <- cbind(
   GT[, c("Sample", "_primary_site")],
   "_primary_disease" = "Normal",
   GT[, c("_gender", "_cohort")]
)
names(GT) <- names(colData)

### Combine, cast factors and order
colData <- rbind(TT, GT)
#### Primary site
colData$`Primary site` <- gsub("Adrenal gland", "Adrenal Gland", colData$`Primary site`)
colData$`Primary site` <- gsub("Cervix Uteri", "Cervix", colData$`Primary site`)
colData$`Primary site` <- gsub("Soft tissue,Bone", "Soft tissue/Bone", colData$`Primary site`)
colData$`Primary site` <- gsub("<not provided>", "", colData$`Primary site`)
colData$`Primary site`[colData$`Primary site` == ""] <- NA
colData$`Primary site` <- as.factor(colData$`Primary site`)
#### Gender
colData$Gender[colData$Gender == ""] <- NA
colData$Gender <- factor(colData$Gender, levels = c("female", "male"))
#### Study
colData$Study <- as.factor(colData$Study)
#### Order and rownames
colData <- colData[order(colData$Study, colData$`Primary site`, colData$`Sample ID`), ]
rownames(colData) <- colData$`Sample ID`



## rowData - genes

### TCGA/TARGET (same as GTEX?)
genes <- read.table("TCGA/gencode.v23.annotation.gene.probemap", sep="\t", header = T, check.names = F, stringsAsFactors = F)


## Read transcript data in, log2(TPM + 0.001)

library(data.table)

### TCGA/TARGET
# TCGA <- fread("TCGA/tcga_target_no_normal_rsem_gene_tpm.gz")
TCGA <- read.table("TCGA/tcga_target_no_normal_rsem_gene_tpm.gz", sep="\t", header = T, check.names = F, stringsAsFactors = F)
rownames(TCGA) <- TCGA$sample
TCGA$sample <- NULL
TCGA <- as.matrix(TCGA)

### GTEX
GTEX <- read.table("GTEX/gtex_RSEM_gene_tpm.gz", sep="\t", header = T, check.names = F, stringsAsFactors = F)
rownames(GTEX) <- GTEX$sample
GTEX$sample <- NULL
GTEX <- as.matrix(GTEX)

### Merge
identical(rownames(TCGA), rownames(GTEX))  # TRUE

TOIL <- cbind(GTEX, TCGA)
# table(colData$`Sample ID` %in% colnames(TOIL))
# table(colnames(TOIL) %in% colData$`Sample ID`)  # only one is missing

colData <- colData[colData$`Sample ID` %in% colnames(TOIL), ]
TOIL <- TOIL[genes$id, colData$`Sample ID`]  # order by genes and colData


### Genes


## Save data as tsv files
dir.create("combined")
write.table(colData, "combined/TOIL-GTEX_TARGET_TCGA-COLUMN_DATA.tsv", sep = "\t", col.names = T, row.names = F, na="", quote=F)
write.table(genes, "combined/TOIL-GTEX_TARGET_TCGA-ROW_DATA.tsv", col.names = T, row.names = F, na="", quote=F)
write.table(TOIL, "combined/TOIL-GTEX_TARGET_TCGA-ASSAY_DATA-LOG2TPM.tsv", col.names = T, row.names = F, na="", quote=F)

# on 03v
# system("gzip -9 TOIL-GTEX_TARGET_TCGA-ASSAY_DATA-LOG2TPM.tsv")
# system("gzip -9 TOIL-GTEX_TARGET_TCGA-COLUMN_DATA.tsv")
# system("gzip -9 TOIL-GTEX_TARGET_TCGA-ROW_DATA.tsv")

##############  Summarized experiment

library(SummarizedExperiment)

rowRanges <- GRanges(
   seqnames = genes$chrom,
   ranges = IRanges(
      start = genes$chromStart,
      end = genes$chromEnd
   ),
   strand = genes$strand,
   feature_id = genes$id
)


TOIL_SE <- SummarizedExperiment(assays=list(TPM=TOIL),
                               rowRanges=rowRanges, colData=colData)

saveRDS(TOIL_SE, "combined/TOIL-GTEX_TARGET_TCGA-SE.Rds")


rm(TCGA, GTEX, GT, TT)
