library(data.table)
source("/home/ev250/Genotyping_RNA_seq/Functions/name_vcf.R") # name txt files made from vcf files
source('/home/ev250/Cincinatti/Functions/various.R')
source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

#' Select variants called with DP >= set value for variants called using RNA
#'
#' @param body full file name to vcf body file
#' @param x cut-off for DP value to call genotypes
#' @param out full name to output file with vcf body to use in shapeit
#' @export
#' @return saves vcf file body to file
#' DP_sel()

DP_sel <- function(body, x, out){
    dt <- name(body)
    gt.col <- grep("_GT", names(dt), value=T)
    dp.col <- grep("_DP", names(dt), value=T)
    ## replace GT based on DP per sample
    for(i in seq_along(gt.col)){
        dt[get(dp.col[i]) < x, gt.col[i] := "./."]
    }
    ## prepare body for vcf to use in shapeit
    ## add ID col:
    dt[, ID :="."]
    vcf <- dt[, c("CHROM","POS","ID", "REF","ALT","QUAL", gt.col),with=F]
    ## add required col FILTER, also phaser defaults to FILTER column in vcf to be PASS
    vcf[,FILTER:="PASS"][,INFO:="."][,FORMAT:="GT"]

    ## reorder cols
    setcolorder(vcf, c("CHROM","POS","ID", "REF","ALT","QUAL", "FILTER", "INFO", "FORMAT", gt.col))
    ## save as tab:
    write.table(vcf, quote=F, row.names=F, file=out, col.names=F, sep = "\t" )

}


## running function
DP_sel(body=snakemake@input[['body']],x=snakemake@params[['DP']], out=snakemake@output[['body']])
