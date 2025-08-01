# Evaluation of antigen pair candidate

# Import libraries
library(data.table)
library(grid)
library(rpart)
library(tools)

# Import data
args <- commandArgs(trailingOnly = TRUE)
indication <- as.character(args[1])
# indication <- "ovarian_serous_cystadenocarcinoma"
indication <- gsub("_", " ", indication)

processed_data_path <- "data/processed/"

# TXT file of processed RNA-seq data 
sheet <- data.table::fread(file=paste0(processed_data_path, "TCGA_GTEX_integrated_antigen_searching_space_selected.txt"))

# TXT file of geometric sketching result
gs <- data.table::fread(file=paste0(processed_data_path, "sketches/", gsub(" ", "-", indication), "-sketch.txt"), header = F)

# CSV file of single antigen candidate list
single.cta <- data.table::fread(file=paste0(processed_data_path, "intermediate/", gsub(" ", "_", indication), "/", "single_ct_antigen.csv"))
single.csfm <- data.table::fread(file=paste0(processed_data_path, "intermediate/", gsub(" ", "_", indication), "/", "single_cellsurface_antigen.csv"))

# CSV file of dual antigen pair candidate list is brought into R
double <- data.table::fread(file=paste0(processed_data_path, "intermediate/", gsub(" ", "_", indication), "/", "dual_antigen.csv"))

# Split the data into test and training 
# target antigen list
antigen.list <- unique(c(sub("(.*?):(.*$)", "\\1", double$combo, perl = T), sub("(.*?):(.*$)", "\\2", double$combo, perl = T), single$combo))

sheet <- sheet[, c("Sample.ID", "type", "tissue.cancer", antigen.list), with=F]
sheet <- sheet[type=="normal"|tissue.cancer==indication]

# split the data set
train <- sheet[Sample.ID%in%gs$V1]
test <- sheet[!Sample.ID%in%gs$V1]

# Function to calculate ADAP cut-off
# relevant max 
relevant_max <- function(y){
    if(IQR(y)==0){
        background <- quantile(y, probs = 0.75)
        y <- y[y>background]
    }
    
    if(length(y)>1){
        return((quantile(y, probs = c(0.75)) + 1.5*IQR(y)))
    }else if(length(y)==1){
        return(y)
    }else{
        return(background)
    }
}

# ADAP cutoff
adap_cutoff <- function(gene){
    # select data
    gene <- as.character(gene)
    tmp <- train[, c("tissue.cancer", gene), with=F]
    tmp <- tmp[tissue.cancer!=indication]
    data.table::setnames(tmp, old=gene, new="gene")
    
    # calculate abs max and rel max for each tissue type
    threshold <- tmp[, list(abs_max=max(gene), rel_max=relevant_max(gene)), by=list(tissue.cancer)]
    
    # calculate cutoff
    cutoff <- min(quantile(threshold$rel_max, probs = 0.95), quantile(threshold$abs_max, probs = 0.95))
    
    return(cutoff)
}

# Use ADAP cutoff line to predict label or tumour vs normal
adap_pred <- function(val, gene){
    gene <- as.character(gene)
    cutoff <- adap_cutoff(gene)
    
    return(data.table::fifelse(val>=cutoff, 1, -1))
}

# Function to find decision tree boundaries on training data
# use the DT cutoff line to predict label of tumour vs normal
predictOnSplits <- function(x, cutoff, ncat)
{
  if (is.null(cutoff))
  {
    return(NA)
  }
  pred = ifelse(x < cutoff, 1, -1)
  if (ncat == -1) # greater than
  {
    pred = ifelse(x > cutoff, 1, -1)
  }
  return (pred)
}

DT_pred <- function(val, gene){
    # select data
    gene <- as.character(gene)
    tmp <- train[, c("type", gene), with=F]
    data.table::setnames(tmp, old=gene, new="gene")
    
    # label
    tmp[, y := fifelse(type=="cancer", 1, -1)]
    tmp[, type := NULL]
    
    # decision tree
    dt = rpart::rpart(y ~ gene, method="class", data=tmp, cp=-1, maxdepth=1)
    
    return(predictOnSplits(x=val, cutoff=dt$splits[4], ncat=dt$splits[2]))
}

# Single antigen predication
x <- data.table::melt(test, id.var=c("Sample.ID", "type", "tissue.cancer"), variable.name="Gene", value.name="Value")

for(i in unique(x$Gene)){
    temp <- x[Gene==i]
    # predict by using adap cut-off on training data
    temp[, Pred_adap := adap_pred(val=Value, gene=i)]
    
    # predict by using decision tree boundaries on the training data
    temp[, Pred_DT := DT_pred(val=Value, gene=i)]
    
    if(exists("res")){
        res <- rbind(res, temp)
    }else{
        res <- temp
    }
}

test.adap <- data.table::dcast(res, Sample.ID+type+tissue.cancer~Gene, value.var="Pred_adap")

test.dt <- data.table::dcast(res, Sample.ID+type+tissue.cancer~Gene, value.var="Pred_DT")

rm(res)

# Function to calculate performance metrics
performance <- function(label, pred){
  tp = sum(label==1&pred==1) # ture positive
  tn = sum(label==(-1)&pred==(-1)) # ture negative
  fp = sum(label==(-1)&pred==1) # false positive
  
  prec = tp / (tp + fp) # precision = tp / (tp + fp)
  rec = tp / sum(label==1) # recall = tp / total number of positives (fn) 
  f1 = (2 * prec * rec) / (prec + rec)
  
  return(list(F1=f1, prec=prec, recall=rec))
}

# Function to evaluate single antigen
evaluationsingle <- function(gene, cutoff){
    gene <- as.character(gene)
    if(cutoff=="ADAP"){
        df <- test.adap[, c("type", gene), with=F]
    }
    
    if(cutoff=="DT"){
        df <- test.dt[, c("type", gene), with=F]
    }
    
    data.table::setnames(df, old=gene, new="gene")
    
    df[, y := fifelse(type=="cancer", 1, -1)]
    
    res <- performance(label = df$y, pred = df$gene)
    
    return(list(Antigen=gene, F1=res$F1, prec = res$prec, recall = res$recall))
}

# Function to evaluate double antigen pairs
evaluationdouble <- function(genes, cutoff){
    genes <- c(as.character(genes["V1"]), as.character(genes["V2"]))
    if(cutoff=="ADAP"){
        df <- test.adap[, c("type", genes[1], genes[2]), with=F]
    }
    
    if(cutoff=="DT"){
        df <- test.dt[, c("type", genes[1], genes[2]), with=F]
    }
    
    data.table::setnames(df, old=c(genes[1], genes[2]), new=c("geneA", "geneB"))
    
    df[, y := fifelse(type=="cancer", 1, -1)]
    
    if(cutoff=="ADAP"){
        df[, Pred := pmax(geneA, geneB)]
    }
    
    if(cutoff=="DT"){
        df[, Pred := pmin(geneA, geneB)]
    }
    
    res <- performance(label = df$y, pred = df$Pred)
    
    return(list(Combo=paste0(genes[1], ":", genes[2]), F1=res$F1, prec = res$prec, recall = res$recall))
}

# Function to calculate the co-expression of two genes (spearman rank)
co_expression <- function(genes){
    genes <- c(as.character(genes["V1"]), as.character(genes["V2"]))
    tmp <- sheet[, c(genes[1], genes[2]), with=F]
    data.table::setnames(tmp, old=c(genes[1], genes[2]), new=c("geneA", "geneB"))
 
    # remove zero expression samples
    zero <- min(tmp$geneA)
    tmp <- tmp[geneA>zero,]
    zero <- min(tmp$geneB)
    tmp <- tmp[geneB>zero,]

    # spearman correlation coefficient
    spearman <- data.table::fifelse(NROW(tmp)>5, cor(tmp$geneA, tmp$geneB, method="spearman"), NaN)

    return(list(Combo=paste0(genes[1], ":", genes[2]), cor=spearman))
}

# Single antigen evaluation
# ct antigen
single.cta[, c( "db", "dist.man", "ttype") := NULL]
# Efficacy
ares <- apply(single.cta, 1, FUN = function (x) {evaluationsingle(gene=x, cutoff = "ADAP") } )
single.cta.adap <- do.call(rbind.data.frame, ares)
single.cta.adap <- data.table::as.data.table(single.cta.adap)
rm(ares)

# tumour specificity
ares <- apply(single.cta, 1, FUN = function (x) {evaluationsingle(gene=x, cutoff = "DT") } )
single.cta.dt <- do.call(rbind.data.frame, ares)
single.cta.dt <- data.table::as.data.table(single.cta.dt)
rm(ares)

# merge two results
data.table::setnames(single.cta.adap, old=c("F1", "prec", "recall"), new=c("OR_F1", "OR_prec", "OR_recall"))
data.table::setnames(single.cta.dt, old=c("F1", "prec", "recall"), new=c("AND_F1", "AND_prec", "AND_recall"))

single.cta.res <- merge(single.cta.adap, single.cta.dt, by="Antigen")
rm(single.cta.adap, single.cta.dt)

# cell surface marker
single.csfm[, c( "db", "dist.man", "ttype") := NULL]
# Efficacy
ares <- apply(single.csfm, 1, FUN = function (x) {evaluationsingle(gene=x, cutoff = "ADAP") } )
single.csfm.adap <- do.call(rbind.data.frame, ares)
single.csfm.adap <- data.table::as.data.table(single.csfm.adap)
rm(ares)

# tumour specificity
ares <- apply(single.csfm, 1, FUN = function (x) {evaluationsingle(gene=x, cutoff = "DT") } )
single.csfm.dt <- do.call(rbind.data.frame, ares)
single.csfm.dt <- data.table::as.data.table(single.cta.dt)
rm(ares)

# merge two results
data.table::setnames(single.csfm.adap, old=c("F1", "prec", "recall"), new=c("OR_F1", "OR_prec", "OR_recall"))
data.table::setnames(single.csfm.dt, old=c("F1", "prec", "recall"), new=c("AND_F1", "AND_prec", "AND_recall"))

single.csfm.res <- merge(single.csfm.adap, single.csfm.dt, by="Antigen")
rm(single.csfm.adap, single.csfm.dt)

# Dual antigen pairs evaluation
double[, V1 := sub("(.*?):(.*$)", "\\1", combo, perl = T)]
double[, V2 := sub("(.*?):(.*$)", "\\2", combo, perl = T)]
double[, c("combo", "db", "dist.man", "ttype") := NULL]

# Efficacy
ares = apply(double, 1, FUN = function(x) {evaluationdouble(genes=x, cutoff = "ADAP")})
double.adap <- do.call(rbind.data.frame, ares)
double.adap = data.table::as.data.table(double.adap)
rm(ares)

# tumour specificity
ares = apply(double, 1, FUN = function(x) {evaluationdouble(genes=x, cutoff = "DT")})
double.dt <- do.call(rbind.data.frame, ares)
double.dt = data.table::as.data.table(double.dt)
rm(ares)

# co-expresssion
ares = apply(double, 1, FUN = function(x) {co_expression(genes=x)})
double.cor <- do.call(rbind.data.frame, ares)
double.cor <- data.table::as.data.table(double.cor)
rm(ares)

# merge two results
data.table::setnames(double.adap, old=c("F1", "prec", "recall"), new=c("OR_F1", "OR_prec", "OR_recall"))
data.table::setnames(double.dt, old=c("F1", "prec", "recall"), new=c("AND_F1", "AND_prec", "AND_recall"))

double.res <- merge(double.adap, double.dt, by="Combo")
double.res <- merge(double.cor, double.res, by="Combo")

# Result output
filename <- paste0("single_antigen_cta_evaluation.csv")
if(!exists(filename)){
    data.table::fwrite(single.cta.res, file=filename)
}

filename <- paste0("single_antigen_CellSurfaceMarker_evaluation.csv")
if(!exists(filename)){
    data.table::fwrite(single.csfm.res, file=filename)
}

filename <- paste0("dual_antigen_evaluation.csv")
if(!exists(filename)){
    data.table::fwrite(double.res, file=filename)
}
