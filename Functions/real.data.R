library(data.table)
source('/home/ev250/Cincinatti/Functions/various.R')


#' Open vcf in R via bcftools query
#'
#' This function allows you to extract part of a vcf file and read it as data table
#' @param cl1 command line bcftools query to extract body of vcf
#' @param cl2 command line bcftools query to extract header of vcf
#' @param sep separation fields based on bcftools query, defaults =" "
#' @keywords vcf fread
#' @export
#' @return data table 
#' vcf_cl()

vcf_cl <- function(cl1,cl2, sep=" ") {
    tmp <- fread(cl1, header=F, sep=sep, colClasses=c(V4="character", V5="character"))
    names(tmp) <- gsub("^.*]","",names(fread(cl2)))[-1]
    return(tmp)
}

#' command line to extract GT and ASE via bcftools query
#'
#' This function allows you to write the command to extract GT and ASE for a region of a vcf
#' @param vcf full path to vcf
#' @param chr chromosome to extract
#' @param st start position to extract
#' @param end end position to extract
#' @param part whether to extract body or header of vcf
#' @keywords command line bcftools query GT ASE
#' @export
#' @return 
#' cl_bcfq()

cl_bcfq<- function(vcf,chr,st, end, part=c("body", "header")) {
    if(part=="body"){
        x <- paste0('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf,  " -r ",chr,":",st,"-",end)
    } else {
        x <- paste0('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf ," -H  | head -1 ")
        }
    
    return(x)
}

#' Extract GT and ASE fro a gene
#'
#' This function allows you to extract GT and ASE for a gene, wraper for cl_bcfq and vcf_cl plus some formatting, removes samples with missing data or unphased. also removes snps if homo or missing for all samples. Adds id col to ease matching with legend/hap format
#' @param vcf full path to vcf
#' @param chr chromosome to extract
#' @param st start position to extract
#' @param end end position to extract
#' @keywords vcf fread
#' @export
#' @return data table 
#' vcf_w()

vcf_w <- function(vcf,chr, st, end) {
    
    body <-cl_bcfq(vcf, chr, st, end , part="body")
   
    header <- cl_bcfq(vcf, chr, st, end , part="header")

    ## open GT and ASE for the selected gene

    gt.ase <- tryCatch({vcf_cl(body,header,sep=" ")}, error=function(e){paste("Region not found", chr,st,end, sep=":")})
    if(length(grep("Region not found", gt.ase))==1){
        return(gt.ase)
    } else {
    
    ## recode names gt.ase to make it compatible with Cincinatti files and functions
    names(gt.ase) <- gsub(":","_",names(gt.ase))

    ## replace unphased data with "." missing value
    gt.ase[gt.ase=="0/0" | gt.ase=="0/1" | gt.ase=="1/0" | gt.ase=="1/1"] <- "."

    ## exclude  snps if homo or missing for all samples
    ex <- apply(gt.ase[,grep("_GT", names(gt.ase)), with=F], 1, function(i) identical(unique(i),c("0|0", ".") ) | identical(unique(i),c("1|1", ".")) |  identical(unique(i),"0|0") | identical(unique(i), ".") | identical(unique(i),"1|1"))

    gt.ase <- gt.ase[which(ex==FALSE),]

    ## this col will help to match snps with legend/hap reference panel
    gt.ase[, id:= paste(POS, REF, ALT, sep=":")]

        return(gt.ase)
    }
    
}




#' get start and end of cis-window per gene to make a bcftools query from a vcf file
#'
#' This function allows you to write the command to extract GT and ASE for a region of a vcf
#' @param file full path to file with gene coordinates, as prepared in inputs.R
#' @param chr chromosome to extract
#' @param gene gene id
#' @param cw length of cis-window, defaults to 5*10^5
#' @keywords command line bcftools query cis window
#' @export
#' @return vector with start and end of cis-window for the selected gene
#' cl_coord()

cl_coord<- function(file,chr,gene,cw=500000){
    ##cat(x) check if command looks ok then run with system, copy cat(x) output to shell and check if it works
    g.st=paste0("awk '$2 ==",chr,"' ",  file, " | grep ", gene, " | cut -d ' ' -f4 | sed 's/,.*//' ")

    g.end=paste0("awk '$2 ==",chr,"' ",  file, " | grep ", gene, " | cut -d ' ' -f5 | sed 's/.*,//' ")

    window.st <- as.numeric(system(g.st, intern=TRUE)) - cw
    
    window.end <- as.numeric(system(g.end, intern=TRUE)) + cw

    v=c(window.st,window.end)
    names(v) <- c("start","end")
    return(v)
}


#' Extract haps from hap file for a subset of snps
#'
#' This function allows you extract haplotypes for a set of snps
#' @param file1 full path to file to legend.gz file
#' @param file2 full path to hap.gz file
#' @param snps vector with POS:REF:ALT coordinates per snps
#' @keywords haplotypes reference panel set snps
#' @export
#' @return matrix with haplotypes as in hap format (each col one hap), rownames are snp id
#' fhaps()

fhaps<- function(file1, file2, snps){
    ## get line number for snp from legend file
    line <- lapply(1:length(snps), function(i) paste0("zcat ", file1, " | sed -n  '/", snps[i], "/='"))  #sed a little bit quicker than awk and grep
    
    nlines <- lapply(line, function(i) as.numeric(system(i, intern=TRUE))) ##  slow
    names(nlines) <- snps
    
    ## sort before reading
    nlines <- sort(unlist(nlines))
    
    ## get haplotypes for fsnps in ref panel, extract from hap.gz file
    ## first line in legend file is headings, need to substract 1 to match hap.gz file
    p.nlines <- paste0(nlines-1, "p")
    
    haps <- paste0("zcat ", file2, " | sed -n '", paste(p.nlines, collapse="; "), "' ")
    rf <- fread(haps, header=F)  ## referene panel for snps in hap format
    mat <- as.matrix(rf)
    rownames(mat) <- names(nlines)
    return(mat)
}


#' Total gene and ASE counts, per fsnp, per individual
#' 
#' Get total and AS counts per snp per individual
#' @param x DT with ASE and GT created from reading vcf
#' @param y data table with total counts for samples
#' @param z data table with each row the genotype for 1 rsnp coded as 0,1,-1 or 2, output from a rec_mytrecase_rSNPs
#' @keywords counts gene ASE 
#' @export
#' @return list of  data tables, each data table corresponds to each rsnp, cols are total counts (y), GT 0,1,-1,2 for the rSNP, ase counts per fsnps across samples
#' tot.ase_counts

tot.ase_counts <- function(x,y=NULL,z=NULL){
    as <- grep("_AS",names(x), value=T)
    tmp <- x[,as,with=F]
    ## get counts for hap2 (n)
    tmp2 <- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], ","), `[[`, 2))))
    ## get counts for hap1+2 (m)
    tmp12<- sapply(1:ncol(tmp), function(i) as.numeric(unlist(lapply(strsplit(tmp[[i]], ","), function(j) sum(as.numeric(j))))))
    if(class(tmp2)=="numeric"){
        tmp3 <- matrix(data=c(tmp2,tmp12), nrow=2, byrow=T)
        rownames(tmp3)=paste0(x$id,c(".n", ".m"))
        colnames(tmp3)=gsub("_AS","",as)
        tmp3 <- t(tmp3)

        } else {
            rownames(tmp2)  <- paste0(x$id,".n")
            rownames(tmp12) <- paste0(x$id, ".m")
            colnames(tmp2) <- colnames(tmp12) <- gsub("_AS","",as)
            tmp3 <- t(rbind(tmp2,tmp12))
            tmp3 <- tmp3[,sort(colnames(tmp3), decreasing=T)]
        }
    if(is.null(z)){
        return(tmp3)
    } else {
    
    ## add gene counts to ase
    tmp4 <- cbind(t(y), tmp3)
    ## add GT of rsnps
    l  <- lapply(1:nrow(z), function(i) {
        tmp5 <- cbind(t(z[i, grep( "_GT", names(z), value=T) , with=F]) , tmp4)
        colnames(tmp5)[1:2] <- c("rsnp", "y")
        rownames(tmp5) <- gsub("_GT","", rownames(tmp4))
        tmp5 <- data.table(tmp5, keep.rownames=T)
        return(tmp5)}
        )      
        return(l)
    }
    
        

}

#' Get the top rsnp for a gene plus number of fsnps
#'
#' This function allows you extract top rsnp for a set of genes and count number of fsnps
#' @param DT1 data table with sig eqtl per gene, got it from Chris Geuvadis /scratch/wallace/Geuvadis_sig_eqtl
#' @param DT2 data table with fsnps for a chr
#' @param chr chr to extract data
#' @keywords top rsnps fsnps gene chr
#' @export
#' @return DT with number of fsnps, top rsnp for genes in a chr
#' top.gene()

top.gene<- function(DT1,DT2,chr ){
    ## get top sig snp per gene with strongest effect size
    setkey(DT1, Gene, Pval, Beta)
    sig <- DT1[CHROM==chr,][,gene_id:=gsub("\\..*", "",Gene)]
    ## get number of fsnps per gene for genes in sig.genes
    fn <- DT2[,.(fsnp.N=.N) ,by=gene_id][gene_id %in% sig$gene_id,]

    ## add max pval and coord of top rsnps to genes with fsnps from chr22
    ## get top sig snp per gene with strongest effect size

    setkey(sig, Gene, Pval, Beta)

    sig.gene.list <- sig[gene_id %in% fn$gene_id, .SD[1], by=gene_id]
    fn.p <- merge(fn, sig.gene.list, by="gene_id")
    setkey(fn.p, Pval, fsnp.N)
    return(fn.p)
}


#' Run linear regression by lm to identify potential cis-qtl
#'
#' This function allows you run lm for log counts vs genotype
#' @param DT1 data table with each row a gene and raw counts for each sample, needs gene_id col, can have extra cols
#' @param DT2 data table with each row a gene,GT for samples and id for rsnp, samples in same order as DT1, gene in col gene_id and snp id in col id. genes do not need to sorted as in DT1.
#' @keywords lm rsnp
#' @export
#' @return DT with regression outcome, gene id and rsnp id.
#' lm.rsnp()

lm.rsnp <- function(DT1,DT2){
    ret <- list()
    for(i in gt.rs$gene_id){
        log.c <- log(unname(unlist(counts.s[gene_id==i,-"gene_id"])))
        g <- unname(unlist(gt.rs[gene_id==i,grep("_GT", names(gt.rs),value=T), with=F]))
    
        fit <- lm(log.c ~ g)

        ret[[i]] <- data.table(summary(fit)$coeff)[-1,][,SNP:=gt.rs[gene_id==i,id]][,gene_id:=i]
    }
    tmp <- rbindlist(ret)
    setkey(tmp,`Pr(>|t|)`)
    return(tmp)
}

#' Convert summary output into data table
#'
#' This function allows you to convert a summary table to data table for input for xtable
#' @param x object (vector) to make summary from
#' @keywords summary data tab;e
#' @export
#' @return DT with summary info
#' sum.fn()
sum.fn= function (x){
  if( !is.numeric(x)) {stop("Supplied X is not numeric")}
  mysummary = data.frame(
            "Min." =as.numeric( min(x)),
            "1st Qu." = quantile(x)[2],
            "Median" = median(x),
            "Mean" = mean(x),
            "3rd Qu." = quantile(x)[4],
            "Max." = max(x),
            row.names=""  
  )
  
  names(mysummary) = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
  return( mysummary )
  }


############ Functions for working with unknown rSNP GT

#' Get haps of fsnp plus any rsnp per sample
#' 
#' extract haps from vcf file for fsnp and add compatible rsnp based on rp
#' @param x data table of fsnps, rows fsnps and columns phased GT per sample, plus any additional
#' @keywords sample haps unknown rsnp
#' @export
#' @return list with 2 elements: hap1 and hap2, each is a  matrix with cols corresponding to hap GT for each fsnp and rows samples.
#' hap_sam.f

hap_sam.f <- function(x){
    ## get haps for fsnps for each sample
    l <- hap_sam(x) ## from Cincinatti/various.R, list of mat1 hap1 and mat2 hap2
    if(nrow(l[[1]])==1){  ## transpose if 1 row to get same format as when many fsnps
        l <- lapply(l,t)
    }
    for(i in seq_along(l)){ rownames(l[[i]]) <- gsub("_GT","", grep("_GT",names(x), value=T)) }
    return(l)
}


    

    
    
    
    
