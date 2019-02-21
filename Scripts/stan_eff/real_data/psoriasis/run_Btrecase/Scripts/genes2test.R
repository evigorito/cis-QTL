library(data.table)
library(xlsx)
library(biomaRt)

#' Select genes to run with Btrecase
#'
#' @param gwas full name to file with gwas hits and their proximal genes
#' @param drg full name to file with differentially regulated genes in psoriasis vs normal skin
#' @param colclass character vector with colclass to use when reading excel file
#' @param d max distance from SNP to gene
#' @param rpkm cut-off for expression on psoriatic skin to consider DRG
#' @param fc fold change cut-off to consider DRG, default >1
#' @param out full path to output file
#' @export
#' @return saves a txt with gene_id and chromosome
#' sel_genes()

sel_genes <- function(gwas, drg, colclass="character",d,rpkm, fc=1, out){
    ## open drg and select genes to use
    drg <- data.table(read.xlsx2(drg, 1, colClasses=colclass))
    drg <- drg[CaseMedian >= rpkm & FC > fc & RankP<=10^-6,]
    ## get gene_id from biomart
    mart <-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
    g.id <- data.table(getBM(attributes=c("ensembl_gene_id","external_gene_name","chromosome_name"),
                                filters="external_gene_name",
                                values=drg[["Gene.Symbol"]],
                             mart=mart))
    
    drg[g.id[,1:2,with=F], gene_id:=i.ensembl_gene_id, on = c(Gene.Symbol="external_gene_name")]
    drg[g.id, (c( "gene_id", "chrom")) := mget(c( "ensembl_gene_id", "chromosome_name")), on = c(Gene.Symbol="external_gene_name")]
    ## alternative syntax
    # drg[g.id, c( "gene_id", "chrom") := list(i.ensembl_gene_id, i.chromosome_name), on = c(Gene.Symbol="external_gene_name")]

    ## select genes in chrom 1-22
    drg <- drg[chrom %in% 1:22,]
    setkey(drg, chrom)
       
    ## gwas side
    gwas <- fread(gwas)
    gwas <- gwas[gene_dist<=d,]
    setkey(gwas,CHROM,POS, gene_dist)
       
    ## combine gwas with drg, save gene_id an CHROM
    tmp <- unique(rbindlist(list(gwas[,.(gene_id, CHROM)], drg[, .(gene_id, chrom)])))
    setkey(tmp,CHROM)
    write.table(tmp, file=out, row.names=F, quote=F)
}
 
sel_genes(snakemake@input[['gwas']],
          snakemake@input[['drg']],
          colclass=snakemake@params[['colclass']],
          d=snakemake@params[['dist']],
          rpkm=snakemake@params[['RPKM']] ,
          out=snakemake@output[[1]])


