library(data.table)
library(GenomicAlignments)
library(GenomicFeatures)
library(biomaRt)

#' Get genes within pre-defined distance to a group of SNPs
#'
#' @param tables character vector with full names to csv with SNPs to search, must have Chrom and position columns
#' @param gene_coord full name to file with gene name , start and end
#' @param d max distance from snp to gene start or end
#' @param out full path to output file
#' @export
#' @return saves a text file with CHROM, POS, rs_id, Alleles, gene_id, gene_name, gene_start,gene_end gene_dist(SNP->gene) and from which input table the data is coming from
#' close_genes()

close_genes <- function(tables, gene_coord, d, out){
    snps <- lapply(tables, fread)
    
    ## select columns with Chrom and position
    cols <- lapply(snps, function(i) sapply(c("chr", "pos"), function(j) grep(j,names(i), ignore.case=T, value=T)))
    snps.sub <- rbindlist(lapply(seq_along(cols), function(i) snps[[i]][ , cols[[i]], with=F]), idcol="table")
    setkey(snps.sub, Chr, Pos)
    ## make GRanges object
    snps_granges=GenomicRanges::GRanges(seqnames = snps.sub$Chr, 
                                          ranges = IRanges::IRanges(start = snps.sub$Pos,
                                                                    end = snps.sub$Pos))
    values(snps_granges) <- DataFrame(table=snps.sub[["table"]])

    ## make GRanges object with gene_coord +- d to start and end
    gene <- fread(gene_coord)
    gene_granges = GenomicRanges::GRanges(seqnames = gene$chrom, 
                                          ranges = IRanges::IRanges(start = gene$start-d,
                                                                    end = gene$end+d))
    values(gene_granges)  <- DataFrame(gene_id = gene[['gene_id']])

    ## find overlaps:
    gene_olaps = GenomicRanges::findOverlaps(gene_granges, snps_granges, ignore.strand=TRUE)
    gene.side <- as.data.table(gene_granges[queryHits(gene_olaps)])
    
    snp.side <- as.data.table(snps_granges[subjectHits(gene_olaps)])

    ## merge gene and snps and format
    comb <- cbind(snp.side[, .(seqnames, start,table)], gene.side[,.(gene_id) ])
    comb[, seqnames := as.integer(as.character(seqnames))]
    setkey(comb, seqnames, start)
    setkey(snps[[1]], Chr, Pos)
    comb[snps[[1]], Marker := i.Marker]

    ## get info from biomaRt
    mart <-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
    g.names <- data.table(getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                                filters="ensembl_gene_id",
                                values=comb$gene_id,
                                mart=mart))
    comb[g.names, gene_name:=i.external_gene_name, on = c(gene_id="ensembl_gene_id")]

    mart2 <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org",dataset="hsapiens_snp")
    ## get rsids and alleles
    s.info <- rbindlist(lapply(unique(comb$start), function(i) {
        data.table(getBM(attributes=c("refsnp_id", "allele", "chr_name","chrom_start"),
                         filters=c("chr_name", "start", "end"),
                         values=list(chr_name=unique(comb[start==i,seqnames]),
                                     start= i,
                                     end=i),                      
                         mart=mart2))
    }))

    ## get alleles for rsids in snps[[1]] in same format

    rsid <- data.table(getBM(attributes=c("refsnp_id", "allele"),
                             filters="snp_filter",
                             values=snps[[1]]$Marker,
                             mart=mart2
                             ))
    comb[s.info, Marker:= i.refsnp_id, on = c(seqnames="chr_name", start="chrom_start")][s.info, Alleles:= i.allele, on= c(seqnames="chr_name", start="chrom_start")]
    comb[rsid, Alleles:=i.allele, on = c(Marker="refsnp_id")]

    ## add distance to gene
    comb <- merge(comb, gene[,.(gene_id,start,end)] , by="gene_id", all.x=T)
    comb[ , g.st:=start.y-start.x][, g.end:=end-start.x]
    comb[ , gene_dist:=g.st]
    ## recode gene_dist , 0 in gene or g.end if SNP is downstream of gene
    comb[(sign(g.st) == 1 & sign(g.end) == -1) | (sign(g.st) == -1 & sign(g.end) == 1) , gene_dist:=0]
    comb[sign(g.st) == -1 & sign(g.end) == -1,  gene_dist:= abs(g.end)]
    ## remove aux cols and rename
    comb[, c("g.st" , "g.end"):=NULL]
    setnames(comb, c("seqnames", "start.x", "start.y", "end", "Marker"), c("CHROM", "POS", "gene_start", "gene_end", "rs_id"))
    setcolorder(comb, c("CHROM", "POS", "rs_id", "Alleles", "gene_id", "gene_name", "gene_start", "gene_end", "gene_dist"))
     
    write.table(comb, file=out, row.names=F)
}
 

close_genes(tables=snakemake@input[['tables']], gene_coord=snakemake@input[['gene_coord']], d=as.numeric(snakemake@params[['cis_window']]), out=snakemake@output[[1]])
