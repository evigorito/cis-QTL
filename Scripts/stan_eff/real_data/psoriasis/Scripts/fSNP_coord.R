library(data.table)
library("GenomicFeatures")
library("GenomicAlignments")

#' Map fSNPs to genes 
#'
#' @param ebg genomics ranges object exon by genes to extract exon coordinates
#' @param fsnp full name to vcf body with SNP chr, pos, ref and alt
#' @param head_fsnp full name to vcf header for fsnp
#' @param out full name to save text file with fSNP coordinates mapped to gene_id
#'
#' @return saves file with mapping info of fSNPs to genes
#' @export
#' map.fsnps()

map.fsnps <- function(ebg, fsnp, head_fsnp, out){
    rebg <- unlist(reduce(readRDS(ebg)))
    fsnp <- name(fsnp)
    fsnp[ , ID:=paste(CHROM, POS, REF, ALT, sep=":")]
    fsnp_granges = GenomicRanges::GRanges(seqnames = fsnp$CHROM, 
                                          ranges = IRanges::IRanges(start = fsnp$POS,
                                                                    end = fsnp$POS))
    values(fsnp_granges)  <- DataFrame(ID = fsnp[['ID']])
    feature_olaps = GenomicRanges::findOverlaps(rebg, fsnp_granges, ignore.strand=TRUE)
    ## extract matches and rename cols
    genes <- as.data.table(rebg[queryHits(feature_olaps)])
    genes[, gene_id:= rebg[queryHits(feature_olaps)]@ranges@NAMES ]
    snps <- as.data.table(fsnp_granges[subjectHits(feature_olaps)])
    ## add ref and alt alleles to snp
    snps <- merge(snps, fsnp, by="ID", all.x=T, sort=F)
    ## select columns
    comb <- cbind(snps[,.(CHROM, POS, ID, REF, ALT)], genes[, .(gene_id)] )
    write.table(comb, file=out, row.names=F)
    
}


#' Format txt converted vcf files: one file per chr
#'
#' This function allows you to get good format from vcf text.
#' @param file file name with full path to vcf
#' @param head path and file with header, defaults to prefix.header.txt
#' @keywords vcf format
#' @export
#' @return named list with info from vcf files
#' name()
name <- function(file, head=NULL){
    
    temp <- fread(file, sep='\t')
    if(is.null(head)){
        temp2 <- fread(paste0(gsub(".txt","",file),".header.txt"))
    } else {
        temp2 <- fread(head, sep="\t")
    }
    
    names(temp)<- gsub("^.*]","",names(temp2))
    names(temp)<-gsub(":","_", names(temp))
    return(temp)        
}


map.fsnps(ebg=snakemake@input[['ebg']], fsnp=snakemake@input[['fsnps']],
          head_fsnp=snakemake@input[['head']], out=snakemake@output[['fsnps']])
