library(data.table)
library(GenomicAlignments)
library(GenomicFeatures)
library(biomaRt)

#' prepares exon by gene input file for calculating counts per gene and prepares data table gene_coordinates, chromosome, GC content, gene length
#'
#' @param gtf gtf file with genome reference annotations
#' @param out.gene.coord full name to save text file with gene coordinates
#' @return saves a  text file with gene_id start and end, chromosome, GC content (%) and length of longest transcript per gene
#' @export
#'
#' gtf_various()
gtf_various<- function(gtf,out.gene.coord){

    txdb <- makeTxDbFromGFF(gtf, format="gtf", circ_seqs=character())
    ebg <- exonsBy(txdb, by="gene")

    gene_id <- names(ebg)
    st <- min(start(ebg))
    end <- max(end(ebg))
    dt <- data.table(gene_id=gene_id, start=st, end=end)
    gene_id <- names(ebg)
    st <- min(start(ebg))
    end <- max(end(ebg))
    dt <- data.table(gene_id=gene_id, start=st, end=end)

    ## add chrom to dt:
    chrom <- select(txdb, keys= gene_id, columns = "TXCHROM", keytype = "GENEID")

    dt <- merge(dt, chrom, by.x="gene_id", by.y="GENEID", sort=F)
    setnames(dt, "TXCHROM", "chrom")
    

    ## add GC and max transcript length
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

    gc <- data.table(getBM(attributes = c("ensembl_gene_id", "transcript_length", "percentage_gene_gc_content"),
                           filters = "ensembl_gene_id",
                           values = dt$gene_id,
                           mart = ensembl))

    setkey(gc, ensembl_gene_id, transcript_length)

    ## get gc and length for the longest transcript
    gc <- gc[, .SD[.N], ensembl_gene_id]
    setnames(gc, names(gc), c("gene_id", "length", "percentage_gc_content"))

    ## merge gc and dt

    tmp <- merge(dt, gc, by="gene_id")

    ## keep only numeric chromosomes and transform to columns to numeric

    ##tmp <- tmp[chrom %in% 1:22,]
    ##tmp[, (names(tmp)[2:ncol(tmp)]) := lapply(.SD, as.numeric), .SDcols=names(tmp)[2:ncol(tmp)]]
    
    write.table(tmp, file=out.gene.coord, row.names=FALSE)
    
}

## running function
gtf_various(snakemake@input[[1]], snakemake@output[[1]])

