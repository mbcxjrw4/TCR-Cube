# Antigen pair identification

# Import libraries
library(data.table)
library(grid)
library(tools)

# Import data
args <- commandArgs(trailingOnly = TRUE)
indication <- as.character(args[1])
# indication <- "ovarian_serous_cystadenocarcinoma"
indication <- gsub("_", " ", indication)

input_data_path <- "data/input/"
processed_data_path <- "data/processed/"

# TXT file of processed RNA-seq data  
sheet <- data.table::fread(file=paste0(processed_data_path, "TCGA_GTEX_integrated_antigen_searching_space_selected.txt"))

# TXT file of geometric sketching result
gs <- data.table::fread(file=paste0(processed_data_path, "sketches/", gsub(" ", "-", indication), "-sketch.txt"), header = F)

# CSV file of CT antigen list
cta <- data.table::fread(file=paste0(input_data_path, "ct_antigens.csv"))

# CSV flie of cell surface marker list 
csfm <- data.table::fread(file=paste0(input_data_path, "membrane_antigens.csv"))
csfm[, Identifier := NULL]
data.table::setnames(csfm, old = "Gene", new = "CellSurfaceMarker")
csfm <- csfm[!(CellSurfaceMarker %in% cta$CTantigen)]

# data selection based on geometric sketching
sheet <- sheet[Sample.ID%in%gs$V1]

## gene has a high expression on indication: 50% quantile > min 
tumour <- sheet[tissue.cancer==indication]
antigen.exp <- sapply(tumour[, 4:NCOL(tumour)], function(x) (return(quantile(x, probs = c(0.50))>min(x))))
antigen.exp <- colnames(tumour)[4:NCOL(tumour)][antigen.exp]

cta <- cta[CTantigen %in% antigen.exp]
csfm <- csfm[CellSurfaceMarker %in% antigen.exp]

feature <- c("Sample.ID", "type", "tissue.cancer", antigen.exp)
sheet <- sheet[, ..feature]

rm(tumour, antigen.exp, feature, gs)

# Function to calculate Davies–Bouldin index
index.DB.fixed <- function(x,cl,d=NULL,centrotypes="centroids",p=2,q=2){
if(sum(c("centroids","medoids")==centrotypes)==0)
  stop("Wrong centrotypes argument")
  if("medoids"==centrotypes && is.null(d))
    stop("For argument centrotypes = 'medoids' d cannot be null")
  if(!is.null(d)){
    if(!is.matrix(d)){
      d<-as.matrix(d)
    }
    row.names(d)<-row.names(x)
  }
  if(is.null(dim(x))){
    dim(x)<-c(length(x),1)
  }
  x<-as.matrix(x)
  n <- length(cl)
  k <- max(cl)
  dAm<-d
  centers<-matrix(nrow=k,ncol=ncol(x))
  if (centrotypes=="centroids"){
    for(i in 1:k){
      for(j in 1:ncol(x)){
        centers[i,j]<-mean(x[cl==i,j])
      }
    }
  }
  else if (centrotypes=="medoids"){
    for (i in 1:k){
      clAi<-dAm[cl==i,cl==i]
      if (is.null(clAi)){
        centers[i,]<-NULL
      }
      else{
        centers[i,]<-.medoid(x[cl==i,],dAm[cl==i,cl==i])
      }
    }   
  }
  else{
    stop("wrong centrotypes argument")
  }
  S<-rep(0,k)
  for(i in 1:k){                             # For every cluster
    ind <- (cl==i)
    if (sum(ind)>1){
      centerI<-centers[i,]
      centerI<-rep(centerI,sum(ind))
      centerI<-matrix(centerI,nrow=sum(ind),ncol=ncol(x),byrow=TRUE)
      S[i] <- mean(apply((x[ind,] - centerI)^2,1,sum)^q)^(1/q)
      
    }
    else
      S[i] <- 0                         
  }
  M<-as.matrix(dist(centers, method="minkowski", p=p))
  R <- array(Inf,c(k,k))
  r = rep(0,k)
  for (i in 1:k){
    for (j in 1:k){
      R[i,j] = (S[i] + S[j])/M[i,j]
    }
    r[i] = max(R[i,][is.finite(R[i,])])
  } 
  DB = mean(r[is.finite(r)])        
  resul<-list(DB=DB,r=r,R=R,d=M,S=S,centers=centers)
  resul
}

# Function to calculate the clustering scores for a single gene

calcscoresingle <- function(gene) {
  gene <- as.character(gene)
  # get cluster values for a single gene
  tmp <- sheet[,c(gene, "tissue.cancer"), with=F]
  tmp[, target :=  fifelse(tissue.cancer == indication, "target", "other")]
  clust <- data.table::fifelse(tmp$target == "target", 1, 2)
  m <- as.matrix(tmp[,1, with=FALSE])
  mode(m) <- "numeric"

  # get mean expression of tumour and normal tissue
  data.table::setnames(tmp, old=gene, new="gene")
  cm <- tmp[, list(MEAN=mean(gene)), by=list(target)]
  cm <- as.matrix(cm$MEAN)
  cent <- dist(cm, method="manhattan")[1]
  
  # get db-index
  if(cent>0.01){
    db <- index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=1)$DB
  }else{
    db <- NaN
  }
  
  # get H/L type
  fl = "H"
  if (as.numeric(cm[1,]) > as.numeric(cm[2,])){
    fl = "L"
  }
  
  return (list(pair=gene, db, cent, fl))
}

# Function to calculate the clustering scores for an antigen pair
calcscorepair <- function(genes){ 
  genes <- c(as.character(genes["V1"]), as.character(genes["V2"]))
  
  # get cluster values per pair
  tmp <- sheet[,c(genes[1], genes[2], "tissue.cancer"), with=F]
  tmp[, target := fifelse(tissue.cancer == indication, "target", "other")]
  clust <- data.table::fifelse(tmp$target == "target", 1, 2)
  m <- as.matrix(tmp[,1:2, with=FALSE])
  mode(m) = "numeric"
  
  # get mean expression of tumour and normal tissue
  data.table::setnames(tmp, old=genes, new=c("geneA", "geneB"))
  cm <- tmp[, list(MEAN_geneA=mean(geneA), MEAN_geneB=mean(geneB)), by=list(target)]
  cm <- as.matrix(cm[, target := NULL])
  cent <- dist(cm, method="manhattan")[1]

  # get db-index
  if(cent>0.01){
    db <- index.DB.fixed(m, as.integer(clust), centrotypes = "centroids", p=1)$DB
  }else{
    db <- NaN
  }
  # get type (HH, HL / LH, LL)
  fl = "H"
  sl = "H"
  if (as.numeric(cm[1,1]) > as.numeric(cm[2,1])){
    fl = "L"
  }
  if (as.numeric(cm[1,2]) > as.numeric(cm[2,2])){
    sl = "L"
  }
  
  return (list(pair=paste(genes[1], genes[2], sep=":"), db, cent, paste0(fl, sl)))
}

### Single antigen
# clustering scores for every gene's power to separate target from all other normal tissue samples\
# db = Davies Bouldin \
# dist.man = Manhattan distance between the cluster centroids \
# ttype = whether gene is on (H) or off (L) in the target
# for cta
ares <- apply(cta, 1, calcscoresingle)
cta.res <- do.call(rbind.data.frame, ares)
cta.res <- data.table::as.data.table(cta.res)
data.table::setnames(cta.res, old=colnames(cta.res), new=c("combo", "db", "dist.man", "ttype"))

# select high expression candidate
cta.res <- cta.res[ttype=="H"]

# select good separation
cta.res <- cta.res[dist.man>1]

# order by DB index
data.table::setorder(cta.res, db)

# select top 100
cta.res <- cta.res[1:100, ]

# for cell surface marker
ares <- apply(csfm, 1, calcscoresingle)
csfm.res <- do.call(rbind.data.frame, ares)
csfm.res <- data.table::as.data.table(csfm.res)
data.table::setnames(csfm.res, old=colnames(csfm.res), new=c("combo", "db", "dist.man", "ttype"))

# select high expression candidate
csfm.res <- csfm.res[ttype=="H"]

# select good separation
csfm.res <- csfm.res[dist.man>1]

# order by DB index
data.table::setorder(csfm.res, db)

# select top 100
csfm.res <- csfm.res[1:100, ]

### Dual antigen pairs
# clustering scores for every possible antigen pair's ability to separate target samples from all other normal tissue samples
# db = Davies Bouldin \
# dist.man = Manhattan distance between the cluster centroids \
# ttype = whether gene is on (H) or off (L) in the target
combos <- expand.grid(cta$CTantigen, csfm$CellSurfaceMarker)
colnames(combos) <- c("V1", "V2")

# for test
# combos <- combos[1:50, ]

rm(ares)
ares = apply(combos, 1, calcscorepair)
double.res <- do.call(rbind.data.frame, ares)
double.res = data.table::as.data.table(double.res)
data.table::setnames(double.res, old=colnames(double.res), new=c("combo", "db", "dist.man", "ttype"))

# select high expression candidate
double.res <- double.res[ttype=="HH"]

# select good separation
double.res <- double.res[dist.man>2]

# order by DB index
data.table::setorder(double.res, db)

# select top 1000
double.res <- double.res[1:10000, ]

## Result output
filename <- paste0(processed_data_path, "intermediate/single_ct_antigen.csv")
if(!exists(filename)){
    data.table::fwrite(cta.res, file=filename)
}

filename <- paste0(processed_data_path, "intermediate/single_cellsurface_antigen.csv")
if(!exists(filename)){
  data.table::fwrite(csfm.res, file=filename)
}

filename <- paste0(processed_data_path, "intermediate/dual_antigen.csv")
if(!exists(filename)){
    data.table::fwrite(double.res, file=filename)
}
