

source('/home/ev250/Cincinatti/Functions/various.R')


#' Get counts per gene from bam files and gtf file
#' 
#' use STAR bam output to get a DT with raw counts per gene
#' @param egb exons by gene, outout from egb
#' @param samp name of dir with STAR bam file
#' @param bam.name name of bam file, same for all samples as output from STAR, defaults to STAR name
#' @param mode="Union", input for summarizeOverlaps, 
#' @param singleEnd=FALSE, input for summarizeOverlaps, defaults pair end
#' @param ignore.strand=TRUE, summarizeOverlaps
#' @keywords missing genotypes
#' @export
#' @return data table with number of missing genotypes per snp per sample
#' miss_dna

counts_sample_sub<- function(ebg,samp,bam.name="Aligned.sortedByCoord.out.bam",mode="Union",singleEnd=FALSE, ignore.strand=TRUE){
   
    filename <- file.path(samp,bam.name)
    bamfiles <- BamFileList(filename, yieldSize=2000000)
    se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode=mode,
                        singleEnd=singleEnd,
                        ignore.strand=ignore.strand)
    temp <- data.table(assays(se)$counts, keep.rownames=T)
    names(temp) <- c("gene_id", samp)
    ##cat("sample",samp)
    return(temp)
}

## snakemake arguments

ebg <- readRDS(snakemake@input[['ebg']])
path <- dirname(snakemake@input[['bam']])
bam.name <- basename(snakemake@input[['bam']])
mode <- snakemake@params[['mode']]
ignore.strand <- as.logical(snakemake@params[['ignore_strand']])
singleEnd <- ignore.strand <- as.logical(snakemake@params[['singleEnd']])
sample <- basename(path)
out <- snakemake@output[[1]]


## call function
counts <- counts_sample_sub(ebg, path, bam.name, mode, singleEnd, ignore.strand)
names(counts)[2:ncol(counts)] <- sample

## save output
write.table(counts, file=out , row.names=F)
                
