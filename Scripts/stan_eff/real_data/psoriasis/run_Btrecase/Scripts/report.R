#' ---
#' title: Select genes to run with Btrecase from psoriasis data
#' author: Elena
#' output:
#'    pdf_document:
#'     toc: true
#'         
#' ---

## Report built by rule report from
## /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Snakefile

## Load R packages and functions:
library(data.table)
library(xlsx)
library(biomaRt)
library(rmarkdown)
source("/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Functions/aux_report.R")


#' ## Sources
#' * We have GWAS SNPs associated with psoriasis and differentially expressed genes between psoriatic and normal skin.
#' * We want to select genes that are likely regulated by SNPs in cis.


#' ## GWAS SNPs
#' For GWAS hits exported from Nat Commun. 2017; 8: 15382
## For each SNP we searched for genes within 500KB 
## Main characteristics

## open file
gwas <- fread(snakemake@input[['gwas']])
setkey(gwas,CHROM,POS, gene_dist)

#' * Number of unique SNPs
length(unique(gwas$rs_id))

#' * Number of unique genes
length(unique(gwas$gene_id))

#' * Distribution of the number of genes proximal to each SNP by gene-SNP distance
d <- c(0, 1 %o% 10^(c(3,4,5)), 5*10^5)
f(dt=gwas,var1="rs_id", var2="gene_id", c("Number.Genes","Number.SNPs"), d)

#' * Same but excluding SNPs within genes
f(dt=gwas[!rs_id %in% gwas[gene_dist==0,rs_id],] ,var1="rs_id",
  var2="gene_id", c("Number.Genes","Number.SNPs"), d[d>0])


#' ## DEG in psoriatic vs healthy skin
#' * Based on J Invest Dermatol. 2014 Jul;134(7):1828-1838. doi: 10.1038/jid.2014.28. Epub 2014 Jan 17.
#' * Only genes with pval <10^-6 are reported, no selection for fold change.
#' * Gene expression is reported in RPKM, low <1, 1<=med<500, high>=500

## open file
drg=data.table(read.xlsx2(snakemake@input[['drg']],
                          sheetIndex=1 ,
                          colClasses=snakemake@params[['colclass']]))

#' * 7238 DEG with Case Median expression >1
#' * 2465 DEG with Case Median expression >1 & FC >1
#' * Distribution of FC by expression levels in up-regulated genes
#' 
f2(dt=drg[ FC>=1,], var1="CaseMedian", var2="FC",
   range1=c(0, 1, 500, max(drg$CaseMedian)))
 
#' * DEG with high expression in psoriatic skin are mostly keratinocytes expressed genes
#'

#' ## DEG proximal to GWAS hits in healthy or psoriatic skin

#' * Look at the expression leves for genes associated with GWAS hits

## 22 out of 38 GWAS hits within genes are DEG, most with low-medium
## expressin levels in psoriatic skin

f3(drg, gwas, var1="Gene.Symbol",var2="CaseMedian",d=d)

summary(drg[Gene.Symbol %in% gwas[['gene_name']] , CaseMedian ])

#' ## Genes to follow up
#' * Genes within 100KB to GWAS hits
#' * Highly expressed upregulated DEG (RPKM >=500)

#' List of genes: 433 
## Up to 100 KB from gene from gwas and high expressed FC>1 in DRG
unique(c(gwas[gene_dist <= 100000, gene_name], drg[CaseMedian>=500 &FC>1, as.character(Gene.Symbol)]))
