#' ---
#' title: QC for psoriasis data
#' author: Elena
#' date: November, 2018
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---


## Load R packages and functions:

library(data.table)
library(biomaRt)
library(tidyr)
library(ggplot2)
library(rmarkdown)
library(parallel)
source('/home/ev250/Cincinatti/Functions/various.R')
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')


#' ##  STAR alignment
dirs <- list.dirs("/mrc-bsu/scratch/ev250/psoriasis/STAR")
files <- unlist(lapply(dirs, list.files, pattern="Log.final.out",
                       full.names=T))
star <- lapply(files, fread, fill=TRUE, sep="|", header=F)
star <- rbindlist(lapply(star, function(i)
    i[c(6,9:10),][, V2:= as.numeric(gsub( "[^0-9.]+", "", V2))]))

starw <- as.data.table(matrix(star$V2,ncol=length(unique(star$V1)),byrow=T))
names(starw) <- star$V1[seq_along(unique(star$V1))]
sum.star <- apply(starw,2,summary)
colnames(sum.star) <- c("Total reads", "Uniq mapped reads", "Uniq map reads (%)")

sum.star



#' ## Gene expression between normal and psoriatic skin
files=list.files("/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts",
                 pattern=".txt", full.names=T)
names(files) <- sapply(files, function(i) gsub(".txt","",basename(i)))
gexp <- lapply(files, fread)

## transform to long format and merge skin type
gexp <- rbindlist(lapply(gexp, function(i) {
    mat <- as.matrix(i[,2:ncol(i)])
    dt <- data.table(gene_id=rep(i[['gene_id']], ncol(mat)),
                     Reads=as.numeric(mat))
    return(dt)
    }), idcol="Skin")
    

## DRG according to Li et al, 2014, 134(7):1828-1838.
DRG <- c("IFNG", "NOS2", "IL6", "IL24", "IL34")

## get ENS id
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
ens  <- getBM(attributes=c("ensembl_gene_id","external_gene_name" ),
              filters="external_gene_name",
              values=DRG,
              mart=ensembl)

## Plot selected genes
gexpSub <- gexp[gene_id %in% ens$ensembl_gene_id,]
gexpSub <- merge(gexpSub, ens, by.x="gene_id", by.y="ensembl_gene_id")
bp <- ggplot(gexpSub, aes(x=Skin, y=log(Reads + 0.1), group=Skin)) + 
    geom_boxplot(aes(fill=Skin)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

bp + facet_grid(. ~ external_gene_name)


#' ## Variants called by RNA-seq
#' * Variants called by RNA with matching POSITION, REF and ALT allele in the reference panel
#'
#' * Before DP filtering

files <- list.files("/mrc-bsu/scratch/ev250/psoriasis/call_vars/RPvar", pattern="chr[0-9]+_Q20_filtRP.txt",  full.names=T)

rna <- lapply(files, name)

## simplify annotations in rna
rna <- rbindlist(mclapply(rna, ann_vcf, mc.cores = parallel::detectCores()))

## get proportion of samples with missing values "./." per snp
cols.gt <- grep("_GT",names(rna))
rna[, miss.p:= apply(rna, 1, function(i)sum(i=="./.")/length(cols.gt))]

## look by annotation
miss.ann <- rna[, summary(miss.p), ANN]
rna.mat <- matrix(miss.ann[['V1']], ncol=6, byrow=T,
                  dimnames=list(unique(miss.ann$ANN),
                                c("Min", "1st Qu", "Median"," Mean","3rd Qu",    "Max")))
rna.mat <- cbind(rna.mat, rna[,.N,ANN][order(match(ANN,rownames(rna.mat))),][['N']])
colnames(rna.mat)[7] <- "N variants"
rna.mat
#' * After applying DP=10 filter per SNP per samples

## apply filter and update
for(i in cols.gt){
    rna[get(names(rna)[i+1]) < 10, names(rna)[i] := "./."]
}
rna[, miss.p:= apply(rna, 1, function(i)sum(i=="./.")/length(cols.gt))]
miss.ann <- rna[, summary(miss.p), ANN]
rna.mat <- matrix(miss.ann[['V1']], ncol=6, byrow=T,
                  dimnames=list(unique(miss.ann$ANN),
                                c("Min", "1st Qu", "Median"," Mean","3rd Qu",    "Max")))
rna.mat <- cbind(rna.mat, rna[miss.p < 1,.N,ANN][order(match(ANN,rownames(rna.mat))),][['N']])
colnames(rna.mat)[7] <- "N variants"
rna.mat
