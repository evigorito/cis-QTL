library(data.table)
library(rasqualTools)

#' Group counts per samples
#' 
#' makes matrix of counts and log library size adjusted by GC contents as in rasqualtools
#' @param gene_info path to file with gene ID, GC content etc output from rule group_gene_counts
#' @param raw.counts character vector with names of files with raw counts per gene
#' @param counts.out full name to file to save dt of counts
#' @param cov.out full name to file to save matrix with GC adjusted log lib size, rownames matrix are gene_id, cols are samples.
#'
#' group_counts()

group_counts <- function(gene_info, raw.counts, counts.out, cov.out){
    ## process inputs
    lcounts <- lapply(raw.counts, fread)
    counts <- Reduce(function(...) merge(...,by="gene_id"), lcounts)

    ## filter counts, remove genes with 0 counts in all samples
    counts <- counts[which(rowSums(counts[, 2:ncol(counts), with=F])> 0),]

    ## prepare arguments for raqualtools
    counts_matrix <- as.matrix(counts[, 2:ncol(counts), with=F])
    rownames(counts_matrix) <- counts$gene_id

    gene_metadata <- fread(gene_info)
    gene_metadata <- gene_metadata[, .(gene_id, percentage_gc_content)]
    gene_metadata[, percentage_gc_content:=as.numeric(percentage_gc_content)]
    
    size_factors = rasqualCalculateSampleOffsets(counts_matrix, gene_metadata, gc_correct = TRUE)

    ## save files
    write.table(counts, counts.out, row.names=F)
    saveRDS(size_factors, cov.out)
}
      

## get args from snakemake
raw.counts <- unlist(snakemake@input[['counts']])
counts.out <- snakemake@output[[1]]
cov.out <-  snakemake@output[[2]]
gene_info <- snakemake@input[['geneInfo']]

group_counts(gene_info=gene_info, raw.counts=raw.counts, counts.out=counts.out, cov.out=cov.out)



