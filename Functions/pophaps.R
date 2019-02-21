library(data.table)
library(Matrix)

#' Selects SNPs from gene and SNPs or indels from cis-window to add up to a user selected number
#
#' @param gene gene id
#' @param gcoord data table with gene_id, start, end and chrom of genes in 1-based coordinates
#' @param ciswindow int withcis-window length
#' @param n number of snps to select
#' @param haps matrix with rownames snp id and cols haps per individual
#' @param leg data table from legend file, with snps to select from, columns: id (snp id matching rownames haps), position, a0 and a1, exta cols allowed
#' @keywords extract snps gene
#' @export
#' @return data table with snpid (id) as in leg file,position, EAF ALL and EAF AMR
#' sel.snps()

sel.snps <- function(gene,gcoord, ciswindow=200000,n=50, haps, leg) {
    ## get gene-snps
    gene.range <- gcoord[gene_id==gene, .(start,end)]
    gsnps <- leg[TYPE== 'Biallelic_SNP' & position>= gene.range$start & position <= gene.range$end ,.(ID,position, ALL, AMR)]
    if(nrow(gsnps) == n){
        return(gsnps)
    }
    
    if(nrow(gsnps) > n) {
        gsnps <- gsnps[sample(1:nrow(gsnps), n), .(ID,position, ALL, AMR)]
        return(gsnps)
    } 
    n2 = n - nrow(gsnps)
    rsnps <- sample(leg[position < gene.range$start - ciswindow | position > gene.range$end + ciswindow, ID], n2)
    all <- leg[ ID %in% c(gsnps$ID, rsnps), .(ID,position, ALL, AMR)]
    return(all)   
}


#' Extract haps for a range of positions from haps file 
#
#' @param file full name to haps file
#' @param positions vector with start and end position to extract
#' @keywords haplotypes snps position
#' @return sparse matrix of haplotyeps,  with rownames snp id, cols samples
#' hap.pos()

hap.pos <- function(file, positions){
    
 ## get snp info for range including line number, from legend file
    snp.i <- fread(cmd=paste("awk '{if ($3 >= ", positions[1], " && $3 <= " ,positions[2], ") {print $0} }'", file))  ## third field in haps file is POS
    snp.i[, id:= paste(V1,V3,V4,V5, sep=":")]
    snp.i[, V1:=id]
    snp.i[, id:=NULL]
    mat <- Matrix(data.matrix(snp.i[, V6:ncol(snp.i)]), sparse=TRUE)
    rownames(mat) <- snp.i$V1
    return(mat)

}
