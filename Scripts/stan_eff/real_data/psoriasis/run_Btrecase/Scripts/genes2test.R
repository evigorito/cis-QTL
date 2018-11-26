library(data.table)
library(xlsx)
library(biomaRt)

#' Select genes to run with Btrecase
#'
#' @param gwas full name to file with gwas hits and their proximal genes
#' @param drg full name to file with differentially regulated genes in psoriasis vs normal skin
#' @param colclass character vector with colclass to use when reading excel file
#' @param d max distance from SNP to gene
#' @param out full path to output file
#' @export
#' @return saves a txt with gene_id and chromosome
#' sel_genes()

sel_genes <- function(gwas, drg, colclass="character",d, out){
    drg <- data.table(read.xlsx2(drg, 1, colClasses=colclass))
    ## get gene_id from biomart
    mart <-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
    g.id <- data.table(getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                                filters="external_gene_name",
                                values=drg[["Gene.Symbol"]],
                             mart=mart))
    drg[g.id, gene_id:=i.ensembl_gene_id, on = c(Gene.Symbol="external_gene_name")]
    ## gwas side
    gwas <- fread(gwas)
    gwas <- gwas[gene_dist<=d,]
    setkey(gwas,CHROM,POS, gene_dist)

    ## save genes, I start with gwas, save gene_id an CHROM
    tmp <- gwas[,.(gene_id, CHROM)]
    write.table(tmp, file=out, row.names=F, quote=F)
}
 
sel_genes(snakemake@input[['gwas']], snakemake@input[['drg']], colclass=snakemake@params[['colclass']], d=snakemake@params[['dist']], out=snakemake@output[[1]])





