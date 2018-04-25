library(data.table)
library(ggplot2)
library(reshape2)
library("GenomicFeatures")
library("GenomicAlignments")
library(biomaRt)
library(parallel)


source("/home/ev250/Genotyping_RNA_seq/Functions/name_vcf.R") # name txt files made from vcf files
source('/home/ev250/Cincinatti/Functions/various.R')
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

#############################################################
## Prepare files with GT and ASE information
## Chr22
chr22 <- DNAvcf4rasqual(path='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE',pattern="less.tab",chr=22)

## After merging into 1 file with mergevcf as explained in geuv.summary.sh check GT field:

vcf.path <- '/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz'
gt.ase=vcf_w(vcf=vcf.path,qc="yes") ## read vcf
vcf.qc <- vcf.gt.qc(gt.ase)

## to save a new vcf.gz (and tabix) excluding snps or samples with wrong genotype format in vcf.gt.qc , use argumentsvcf.path,  exclude=c("snps","samples") with path (defaults to wd).

exclude.samples <- vcf.gt.qc(gt.ase, exclude="samples",vcf.path, path="/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/", vcf.out="chr22.GTqc")


#################################################################
## Prepare matrix of counts per gene: adapted from ../Cincinatti/Scripts/raw_counts_rnaseq_star37_GA.R

## prepare exon by gene input file for calculating counts per gene 

ebg <- gtf2ebg('/scratch/ev250/reference_genome/built37/Homo_sapiens.GRCh37.87.gtf')

saveRDS(ebg, "/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_ebg.rds")

## run countxgene.R from call_array.sh

## open files and merge

f <- list.files(path="/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/counts_sample",full.names=T)

lcounts <- lapply(f,fread)

counts <- Reduce(function(...) merge(...,by="gene_id"), lcounts)

write.table(counts, '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_raw_counts.txt', row.names=F)

counts <- fread(file='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_raw_counts.txt')

## filter counts ,  mean >= 100 across samples

f <- 100 ## filter
counts.f <- counts[which(rowMeans(counts[, 2:ncol(counts), with=F])>=f),]

write.table(counts.f, '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt', row.names=F)

## if samples were excluded due to QC issues in previous section, those samples need to be excluded from count matrix to run btrecase:

counts.f.ex <- counts.f[,!exclude.samples$samples.excluded, with=FALSE]

write.table(counts.f.ex, '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt', row.names=F)

##########################################################################################################
## Calculate library size: save as matrix in log scale

lib_size= log(colSums(as.matrix(counts[,2:ncol(counts),with=F])))
##save as matrix
lib_size <- matrix(lib_size, ncol=1, dimnames=list(names(counts)[2:ncol(counts)], "lib.s"))
saveRDS(lib_size, '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.rds')

lib_size.GT <- log(colSums(as.matrix(counts[,names(counts.f.ex)[2:ncol(counts.f.ex)],with=F])))
lib.size.GT <- matrix(lib_size.GT, ncol=1, dimnames=list(names(counts.f.ex)[2:ncol(counts.f.ex)], "lib.s"))

saveRDS(lib.size.GT, '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds')

## calculate prior for library size: prepare matrix with cols: ind,counts(per gene, all genes),lib.size(entred and standardised).

lib.sz.sc <- scale(lib_size)
## counts.f
lib.sz.sc.ind <- lib.sz.sc[rep(1:nrow(lib.sz.sc), rep(nrow(counts.f), nrow(lib.sz.sc))),]

counts.long <- melt(counts.f, id.vars="gene_id", measure.vars=names(counts)[2:ncol(counts)])
counts.long[,log.lib.size:=lib.sz.sc.ind]
counts.long[, log.counts:=log(value)]

## linear model: not able to run with more than 1000 samples

fit <- lm(log.counts ~ log.lib.size, data=counts.long)$summary 
##Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
##  NA/NaN/Inf in 'y'
rows <- sample(1:nrow(counts.long),1000)
plot(y=counts.long$log.counts[rows], x=counts.long$log.lib.size[rows])
fit <- summary(lm(log.counts ~ log.lib.size, data=counts.long[rows,]))$coefficients


#################################################################
## select fsnps per gene: get exon info for each gene from biomart

## 1) get exon bounderies for each gene, take the longest exon only
fit
mart <-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

exons <- data.table(getBM(attributes=c("ensembl_gene_id", "chromosome_name","strand","exon_chrom_start","exon_chrom_end"), filters=c("ensembl_gene_id"), values=counts$gene_id, mart=mart))

## make format as required by rasqualtools:

st <- exons[, .(exon_starts=paste0(exon_chrom_start, collapse=",")), by=ensembl_gene_id]
end <- exons[, .(exon_ends=paste0(exon_chrom_end, collapse=",")), by=ensembl_gene_id]

exons.st.end <- merge(st,end, by="ensembl_gene_id")

exons.chr.strand <- exons[,.SD[1], by="ensembl_gene_id"]

exons_all <- merge(exons.chr.strand[,1:3, with=F], exons.st.end, by= "ensembl_gene_id")

#renaming as rasqualTools to avoid errors

names(exons_all) <- c("gene_id","chr","strand", "exon_starts","exon_ends")

## selecting longest exon to avoid double counting
exons_longest <- format_exons(exons_all)

write.table(exons_longest,file="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt", row.names=FALSE, quote=FALSE)

exons_longest <- fread(file="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt")

## 2) get SNP coordinates to define fSNPs and rSNPs: body of vcf file with all samples for chr22

snp.coord <- name('/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.snp.coord.txt', head="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.header.snp.coord.txt")


## 3) get fSNPs per gene

e.22 <- exons_longest[chr==22,]
fsnps.22 <- lapply(1:nrow(e.22), function(i) get.f.Snps2(e.22[i,],snp.coord))
names(fsnps.22) <- e.22$gene_id

## save

saveRDS(fsnps.22, file='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.rds')

## save as txt

write.table(rbindlist(fsnps.22), file='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt', row.names=F)

