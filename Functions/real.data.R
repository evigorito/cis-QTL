
library(data.table)
source('/home/ev250/Cincinatti/Functions/various.R')
##library(snpStats)

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

## #' command line to extract GT and ASE via bcftools query
## #'
## #' This function allows you to write the command to extract GT and ASE for a region of a vcf
## #' @param vcf full path to vcf
## #' @param chr chromosome to extract
## #' @param st start position to extract
## #' @param end end position to extract
## #' @param part whether to extract body or header of vcf
## #' @keywords command line bcftools query GT ASE
## #' @export
## #' @return character vector with bcftools command
## #' cl_bcfq()

## cl_bcfq<- function(vcf,chr=NULL,st=NULL, end=NULL, part=c("body", "header")) {
  
 
##     if(part=="body"){
##         x <- paste0('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf)
##         if(!(is.null(chr) & is.null(st) & is.null(end))){
            
##             x <- paste0('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf,  " -r ",chr,":",st,"-",end)
##         }
        
##     } else {
##         x <- paste0('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf ," -H  | head -1 ")
##         }
    
##     return(x)
## }

#' command line to extract any field via bcftools query, allow to select samples or region if required
#'
#' This function allows you to write the command to extract any field for a region of a vcf, allows sample selection
#' @param vcf full path to vcf or bcf, compressed or uncompressed
#' @param chr chromosome to extract
#' @param st start position to extract
#' @param end end position to extract
#' @param samples character vector with samples to extract, defaults=NULL
#' @param f.arg character vector with -f argument for bcftools query, defaults to GT and ASE
#' @param part whether to extract body or header of vcf
#' @keywords command line bcftools query GT ASE
#' @export
#' @return character vector with bcftools command
#' cl_bcfq()

cl_bcfq<- function(vcf,chr=NULL,st=NULL, end=NULL, samples=NULL, f.arg=NULL, part=c("body", "header")) {

    if(!is.null(samples)){## select samples first
        samp <- paste(" -s ", paste0(samples, collapse=","))
    }
    if(!is.null(chr) & !is.null(st) & !is.null(end)){
        reg <- paste0(" -r ",chr,":",st,"-",end)
    }
    
    
    if(part=="body"){
        if(is.null(samples) & !exists("reg") & is.null(f.arg)){
            x <- paste('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf)
        }
        if(!is.null(samples) & !exists("reg") & is.null(f.arg)){
            x <- paste(' bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', samp, vcf)
        }
        if(is.null(samples) & exists("reg") & is.null(f.arg)){
            x <- paste('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', reg, vcf)
        }
        if(!is.null(samples) & exists("reg") & is.null(f.arg)){
            x <- paste(' bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', reg, samp, vcf)
        }
        if(is.null(samples) & !exists("reg") & !is.null(f.arg)){
            x <- paste('bcftools query -f ', f.arg , vcf)
        }
        
        if(!is.null(samples) & !exists("reg") & !is.null(f.arg)){
            x <- paste(' bcftools query -f ', f.arg, samp, vcf)
        }

        if(!is.null(samples) & exists("reg") & !is.null(f.arg)){
            x <- paste( ' bcftools query -f ', f.arg, reg, samp, vcf)
        }
      
    } else {

        if(is.null(samples) & is.null(f.arg)){
            x <- paste('bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n" ', vcf ," -H  | head -1 ")
        }
        if(!is.null(samples) & is.null(f.arg)){
            x <-  paste0(' bcftools query -f "%CHROM %POS %ID %REF %ALT[ %GT] [ %AS]\\n"',samp, vcf,' -H | head -1')
        }
        if(is.null(samples) & !is.null(f.arg)){
            x <- paste('bcftools query -f  ',f.arg, vcf ," -H  | head -1 ")
        }
        if(!is.null(samples) & !is.null(f.arg)){
            x <-  paste('bcftools query -f', f.arg, samp, vcf,' -H | head -1')
        }
        
    }
    
    return(x)
}




#' Extract GT and ASE for a gene
#'
#' This function allows you to extract GT and ASE for a gene, wraper for cl_bcfq and vcf_cl plus some formatting, removes samples with missing data or unphased. also removes snps if homo or missing for all samples. Adds id col to ease matching with legend/hap format
#' @param vcf full path to vcf
#' @param chr chromosome to extract, null for whole vcf
#' @param st start position to extract, null for whole vcf
#' @param end end position to extract, null for whole vcf
#' @param samples character vector with samples to extract, defaults=NULL
#' @param f.arg character vector with -f argument for bcftools query, defaults to GT and ASE
#' @param qc use function for qc purpose only
#' @param exclude whether to return a list with snps excluded from vcf (all homo or missing)
#' @keywords vcf fread
#' @export
#' @return data table with GT and ASE, unless !is.null(excluded), returns list with first element data table with GT and ASE and second element a data table with exluded snps.
#' vcf_w()

vcf_w <- function(vcf,chr=NULL, st=NULL, end=NULL, samples=NULL, f.arg=NULL, qc=NULL, exclude=NULL) {
    
        body <-cl_bcfq(vcf, chr, st, end ,samples,f.arg, part="body")
   
        header <- cl_bcfq(vcf, chr, st, end , samples,f.arg, part="header")
        

    ## open GT and ASE for the selected gene

    gt.ase <- tryCatch({vcf_cl(body,header,sep=" ")}, error=function(e){paste("Region not found for chrom and positions", chr,st,end, sep=":")})
    if(is.character(gt.ase)){
        return(gt.ase)
    } else {
        
        if(!is.null(qc)){
            return(gt.ase)
            } else {
            
                ## recode names gt.ase to make it compatible with Cincinatti files and functions
                names(gt.ase) <- gsub(":","_",names(gt.ase))

                ## replace unphased data with "." missing value
                gt.ase[gt.ase=="0/0" | gt.ase=="0/1" | gt.ase=="1/0" | gt.ase=="1/1"] <- "."

                ## exclude  snps if homo or missing for all samples
                ex <- apply(gt.ase[,grep("_GT", names(gt.ase)), with=F], 1, function(i) identical(unique(i),c("0|0", ".") ) | identical(unique(i),c("1|1", ".")) |  identical(unique(i),"0|0") | identical(unique(i), ".") | identical(unique(i),"1|1"))

                ## this col will help to match snps with legend/hap reference panel
                gt.ase[, id:= paste(POS, REF, ALT, sep=":")]
                
                if(is.null(exclude)){
                    ## select relevant snps
                    gt.ase <- gt.ase[which(ex==FALSE),]
                    return(gt.ase)
                } else {
                    excl=gt.ase[which(ex==TRUE),]
                    excl[,reason:="Missing or homo GT all samples"]
                    excl <- excl[,.(id,reason)]
                    l <- list(keep= gt.ase[which(ex==FALSE),],excluded=excl)
                    return(l)
                }
        
            }
    
    }
}


#' Check if gt.ase, output from vcf_w has phased GT in GT field, can also save a new vcf with excluding wrong GT by snp or by sample, as required.
#'
#' This function allows you to test wether a vcf_w returned object is correctly formatted in GT field with the option to list snps or samples with wrong GT in format to be excluded from 
#' @param gt.ase object returned from vcf_w
#' @param exclude optional argument,removes entries with wrong GT format by snps or by sample, options "snps" or "samples"
#' @param vcf.path path and file name of original vcf, argument used for preparing new vcf with wrong GT entries removed
#' @param path optional, path to save new vcf excluding wrongly formatted GT, the default corresponds to the working directory
#' @param vcf.out optional, prefix for new vcf with wrong GT entries removed
#' @keywords vcf gt qc 
#' @export
#' @return named vector when the only argument used is gt.ase. The vector gives the total number of snps, number of snps with wrong GT format, total number of samples and number of samples with wrong GT format. When all arguments are used it saves and indexes a new vcf excluding wrongly GT entries by snps or samples in format vcf.gz. In this mode the function returns a DT with the chr:pos:ref:alt of snps excluded or the name of the samples excluded.
#' vcf.gt.qc()

vcf.gt.qc <- function(gt.ase, exclude=c("snps","samples"), vcf.path, path=".", vcf.out="chr22.GTqc") {
 
    gt.col <- grep("GT$", names(gt.ase))
    ok <- c("0|1", "0|0","1|1", "1|0")
    by.snp <- apply(gt.ase[,gt.col,with=FALSE],1, function(i) sum(!i %in% ok)!=0)
    by.sample <- apply(gt.ase[,gt.col,with=FALSE],2, function(i) sum(!i %in% ok)!=0)
    report <- c(nrow(gt.ase), sum(by.snp),length(gt.col), sum(by.sample))
    names(report) <- c("total snps", "snps with wrong GT", "total samples","samples with wrong GT")
    if(!(missing(exclude)  & missing(vcf.path) & missing(vcf.out))) {
        na.ex <- pmatch(exclude,c("snps","samples"))
        if(is.na(na.ex) | missing(vcf.path) | missing(vcf.out)){
            stop("invalid 'exclude','vcf.path' or 'vcf.out' argument")
        } else {
            out <- paste0(path,"/",vcf.out,".vcf.gz")
            if (exclude=="snps"){              
                DT <- gt.ase[which(by.snp),.(CHROM,POS,REF,ALT)]
                DT[,ex:=paste0(CHROM,":",POS)]
                DT[,snps.excluded:=paste0(ex,":",REF,":",ALT)]
                del <- paste0("^",paste0(DT$ex, collapse=","))
                bcf.f <- paste("bcftools view -Oz -t",del,vcf.path,"-o",out)         
                report <- data.table(snps.excluded=DT$snps.excluded)
                
            } else {
                del <- names(by.sample)[by.sample]
                del <- gsub(".GT$","",del)
                del2 <- paste0("^",paste0(del,collapse=","))                
                bcf.f <- paste("bcftools view -Oz -s",del2,vcf.path,"-o",out)
                report <- data.table(samples.excluded=del)
            }
            system(bcf.f)
            bcf.i <- paste("bcftools index -t", out)
            system(bcf.i)            
        }
    }
        
    return(report)
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

    snps <- unique(snps) ## remove duplicate entries
    ## get line number for snp from legend file
    
    l.snps <- paste0("zcat ", file1, " | sed -n '", paste0("/",snps, "/{=;p}", collapse=";"), "'")# | sed '{N;s/\n/ /}'") ## line number first line, snp info second line to keep track of snps if some arent on ref panel
    
    l.s <- system(l.snps, intern=TRUE)

    n.lines <- unlist(as.numeric(l.s[1:length(l.s) %%2 ==1]))
    names(n.lines) <- unname(sapply(sapply(strsplit(l.s[1:length(l.s) %%2 ==0], split=" ") , `[[`,1) , function(i) sub(".*?:(.+)","\\1",i)))

    if(length(n.lines)==0){
        return("no snps in reference panel")

    } else {
    
    
        ## get haplotypes for fsnps in ref panel, extract from hap.gz file
        ## first line in legend file is headings, need to substract 1 to match hap.gz file
        
        p.nlines <- paste0(n.lines-1, "p")
        s.nlines <- n.lines-1

        ## find stretches of sequencial lines to improve running time sed -n '1,4p'
        ## create stretches of 1's
        st <- s.nlines -c(s.nlines[1],s.nlines[-length(s.nlines)])

        if(sum(st==1)>=1){
            v <- rle(unname(st)) ## rle returns an object with two components, lengths and values.         
            w <- which(v$value==1) ## relevant entries from v$value
            vl <- v$length ## to ease coding, I need to sum over elements of vl to get indexes for p.nlines/s.nlines
            
            if(w[1]==2) { ## 2 corresponds to the range of ones starting from the beginning, first element always 0
                com <- paste(s.nlines[w[1]-1],p.nlines[sum(vl[1:w[1]])],  sep=",")
            } else {
                com <- c(paste(p.nlines[1:(sum(vl[1:(w[1]-1)])-1)], collapse="; "), paste(s.nlines[sum(vl[1:(w[1]-1)])],p.nlines[sum(vl[1:w[1]])], sep=","))
                
            }
            if(length(w)>1){       
                for(i in 2:length(w)){
                    if(sum(vl[1:w[i-1]])+vl[w[i]] +1 == sum(vl[1:w[i]])){## the new stretch of 1's starts inmediately after the last stretch
                        com <- c(com, paste(s.nlines[sum(vl[1:(w[i-1]+1)])],p.nlines[sum(vl[1:w[i]])], sep=","))
                    }
                    
                    else if(sum(vl[(w[i-1]+1):w[i]])==2){ ## the new stretch of 1's in w[i] starts one snp after the last stretch of ones
                        com <- c(com,paste(s.nlines[sum(vl[1:(w[i]-1)])],p.nlines[sum(vl[1:w[i]])], sep=","))
                    } else {
                        
                        com <- c(com, c(paste(p.nlines[(sum(vl[1:w[i-1]])+1):(sum(vl[1:w[i]-1])-1)], collapse="; "), paste(s.nlines[sum(vl[1:(w[i]-1)])],p.nlines[sum(vl[1:w[i]])], sep=",")))
                    }
                    
                }             
            }
            ## I may need to add elements after last entry of w
            if(sum(vl[1:w[length(w)]]) < length(p.nlines)) {           
                i=length(w)
                com <- c(com,paste(p.nlines[(sum(vl[1:w[i]])+1):sum(vl)], collapse="; "))             
            }

            haps <- paste0("zcat ", file2, " | sed -n '", paste(com, collapse= "; "), "' ")
            
        } else {              
            
            haps <- paste0("zcat ", file2, " | sed -n '", paste(p.nlines, collapse="; "), "' ") 

        }
        
        rf <- fread(haps, header=F)  ## referene panel for snps in hap format
        mat <- as.matrix(rf)
        rownames(mat) <- names(n.lines)
        return(mat)
    }
}
        
#' Extract haps from hap file for the whole cis-window (range)
#'
#' This function allows you extract haplotypes for a range of snps
#' @param file1 full path to file to legend.gz file
#' @param file2 full path to hap.gz file
#' @param cw vector start and end position of snps to extract within the range
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL
#' @param maf cut-off for maf
#' @keywords haplotypes reference panel range snps
#' @export
#' @return matrix with haplotypes as in hap format (each col one hap), rownames are snp id
#' haps.range()

haps.range<- function(file1, file2, cw,population="EUR", maf=0.05){

    eth <- 6:11 ## field number in legend file for ethniticy
    names(eth) <- c("AFR", "AMR", "EAS",  "EUR", "SAS", "ALL")
    
    ## get snp info for range including line number, from legend file
    snp.i <- system(paste0("zcat ", file1, " | awk '{if ($2 >= ", cw[1], "&& $2 <= " ,cw[2], ") {print NR \" \" $2\":\"$3 \":\" $4 \" \" $",unname(eth[names(eth)==population]),"} }'"), intern=TRUE)  ## second field in legend file is POS,then ref then alt allele

    if(length(snp.i)==0) return("no snps in reference panel")
    
    ## format snp.i as DT
    
    DT <- as.data.table(lapply(1:3, function(i) unlist(lapply(strsplit(snp.i, split=" ") , `[[`,i))))
    names(DT) <- c("line","snp","maf")
    DT[,line:=as.numeric(line)-1]  ## first line in legend file is headings, need to substract 1 to match hap.gz file
    DT[,maf:=as.numeric(maf)]
    haps <- paste0("zcat ", file2, " | sed -n '", DT$line[1], ",", DT$line[nrow(DT)], "p' ")
    
    rf <- fread(haps, header=F)  ## referene panel for snps in hap format
    ## remove snps below maf cut-off
    
    keep <- which(DT$maf>=maf &  DT$maf< (1-maf))
    rf <- rf[keep,]
    mat <- as.matrix(rf)
    rownames(mat) <- DT$snp[keep]
    return(mat)
    
}

#' Extracts snp eaf from legend file for a set of snps
#'
#' This function allows you extract snp info from legend file for a range of snps
#' @param file1 full path to file to legend.gz file
#' @param snps character vector with id of snps to extract information, id=POS:REF:ALT
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL
#' @keywords snp eaf reference panel
#' @export
#' @return DT with snp id and EAF from legend file for selected population
#' snp.eaf()

snp.eaf<- function(file1, snps,population="EUR"){

    eth <- 6:11 ## field number in legend file for ethniticy
    names(eth) <- c("AFR", "AMR", "EAS",  "EUR", "SAS", "ALL")
    ## get first and last POS within snps
    tmp <- as.numeric(gsub(":.*","",snps))
    cw <- c(min(tmp),max(tmp))
    
    ## get snp info for range including line number, from legend file
    snp.i <- system(paste0("zcat ", file1, " | awk '{if ($2 >= ", cw[1], "&& $2 <= " ,cw[2], ") {print NR \" \" $2\":\"$3 \":\" $4 \" \" $",unname(eth[names(eth)==population]),"} }'"), intern=TRUE)  ## second field in legend file is POS,then ref then alt allele

    if(length(snp.i)==0) return("no snps in reference panel")
    
    ## format snp.i as DT
    
    DT <- as.data.table(lapply(1:3, function(i) unlist(lapply(strsplit(snp.i, split=" ") , `[[`,i))))
    names(DT) <- c("line","snp","eaf")
    DT[,eaf:=as.numeric(eaf)]
    DT[,line:=NULL]
    DT <- DT[snp %in% snps,]
    ## DT has unique entries, tmp may have duplicated snps
    DT2 <- merge(data.table(snp=snps),DT, by="snp", sort=F)   
    return(DT2)
}



#' Total gene and ASE counts, per fsnp, per individual
#' 
#' Get total and AS counts per snp per individual
#' @param x DT with ASE and GT created from reading vcf
#' @param y data table with total counts for samples
#' @param z data table with each row the genotype for 1 rsnp coded as 0,1,-1 or 2, output from a rec_mytrecase_rSNPs
#' @keywords counts gene ASE 
#' @export
#' @return matrix if z=NULL or list of  data tables, each data table corresponds to each rsnp, cols are total counts (y), GT 0,1,-1,2 for the rSNP, ase counts per fsnps across samples
#' tot.ase_counts

tot.ase_counts <- function(x,y=NULL,z=NULL){
    as <- grep("_AS",names(x), value=T)
    tmp <- x[,as,with=F]
    ## missing values in phaser for AS are ".", happens when dealing with rna GT, convert to 0,0, as they wont affect counts
    tmp[tmp=="."] <- "0,0"
    
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
        y2 <- y[,which(names(y) %in% rownames(tmp3)), with=FALSE]
        tmp4 <- cbind(t(y2), tmp3)
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

#' Recode GT from trecase scale (0,1,-1,2) to GUESSFM scale 0 M, 1 hom ref, 2 het, 3 hom alt
#'
#' This function allows you to recode GT for input in tag function from GUESSFM
#' @param DT1 data table GT coded in trecase scale, rows SNPS, cols samples plus additionals
#' @keywords recode GUESSFM
#' @export
#' @return Matrix with rows samples and cols SNPS
#' rec.guess()

rec.guess<- function(DT){
    M <- t(as.matrix(DT[, grep("_GT", names(DT),value=T), with=F]))
    ##recode
    M[M==2] <- 3
    M[abs(M)==1] <- 2
    M[M==0] <- 1
    M[is.na(M)] <- 0
    M[M=="."] <- 0
    colnames(M) <- DT$id
    rownames(M) <- grep("_GT", names(DT),value=T)
    M <- apply(M,2,as.numeric)
    return(M)
    
}


#' function for filtering input for stan, also option to remove fsnps with 0 m counts in all samples
#'
#' This function allows you to test whether a rsnp has enough het inds with a certain number of ASE counts
#' @param geno.exp data table for an element of output list tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5
#' @param rem whether to exclude fsnps with 0 m counts in all samples, defaults to NULL
#' @keywords filtering stan input
#' @export
#' @return data table with sample name, GT rsnp, total counts, n and m counts for each fSNP, same data table outputted as an element of tot.ase_counts
#'
#' filt.rsnp()

filt.rsnp <- function(geno.exp,ase=5, n=5, rem=NULL){
   
    n.col <- grep("\\.n",names(geno.exp), value=T)
    m.col <- grep("\\.m",names(geno.exp), value=T)
    m.counts <-  rowSums(geno.exp[,m.col,with=F])
    # select ASE input when total ase counts are above threshold and 
    A <- which(m.counts>=ase)
    if(nrow(geno.exp[A,][abs(rsnp)==1,])<n)
        return("Not enough individuals with ASE counts")
    if(is.null(rem)) return(geno.exp)
    ## remove
    tmp=geno.exp[,m.col, with=FALSE]
    rem=names(tmp)[colSums(geno.exp[,m.col,with=F])==0]
    s=sapply(rem, function(i) sub("\\.m","",i))
    r=sapply(s, function(i) grep(i, names(geno.exp)))
    return(geno.exp[,r:=NULL])
}


#' function for filtering fsnps, allows to select filtering if n=0 or n=m (likely GT error), also option to remove fsnps with less than cut-off total ASE counts in all samples, INPUT matrix
#'
#' This function allows you to test whether a rsnp has enough het inds with a certain number of ASE counts
#' @param c.ase matrix output of tot.ase_counts
#' @param ase total ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.snp per snp ase cut-off
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5
#' @param gt.err optional parameter, if "yes" it sets to m=0 when corresponding n=0 or n=m (likely GT error, homo called het), defaults to NULL
#' @param rem whether to exclude fsnps with less than min.ase.snps ASE counts in all samples, defaults to "yes", otherwise set it to null
#' @keywords filtering matrix m counts
#' @export
#' @return matrix with rownames sample name, cols: n and m counts for each fSNP, option to remove snps with ASE counts below cut-off in all individuals. Also ASE count below threshold (ase.min.snp) are replaced with 0, so will be excluded in subsequent analysis
#'
#' filt.fsnp()

filt.fsnp <- function(c.ase,ase=5, min.ase.snp=5, n=5, gt.err=NULL, rem="yes"){
  
    n.col <- grep("\\.n",colnames(c.ase), value=T)
    m.col <- grep("\\.m",colnames(c.ase), value=T)
    ## remove all entries below cut-off, both for n and m cols
    c.ase[c.ase < min.ase.snp] <- 0  ## entries with less than cut-off converted to 0
    ## remove counts when n=0 or n=m
    if(!is.null(gt.err)){
        for(x in 1:length(n.col)){
        w <- c.ase[,n.col[x]]==0 | c.ase[,n.col[x]]==c.ase[,m.col[x]] 
        c.ase[w,m.col[x]] <- 0
        }
    }
    
    
    m.counts <-  rowSums(c.ase[,m.col, drop=F])
    ## select ASE input when total ase counts are above threshold, only inlcudes snps with counts > min.ase.snps
    A <- which(m.counts>=ase)
    if(nrow(c.ase[A,,drop=FALSE]) < n ) return("Not enough individuals with sufficient ASE counts per exonic snp")
    ## remove
    if(is.null(rem)) return(c.ase)
    tmp = c.ase[,m.col]
    keep=colnames(tmp)[colSums(tmp)!=0]  ## 
    s=sapply(keep, function(i) sub("\\.m","",i))
    k= sapply(s, function(i) grep(i, colnames(c.ase), value=T))
    return(c.ase[, k])
}



#' quick association test  for filtering out rnsps as input for stan
#'
#' This function allows you to predict whether an rsnp is likely to be associated with gene expression
#' @param geno.exp data table for an element of output list tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trecase default 5
#' @keywords filtering stan input
#' @export
#' @return vector with correlation of y and g and z.test for p
#'
#' q.rsnp()

q.rsnp <- function(geno.exp,ase=5, n=5){
   
    n.col <- grep("\\.n",names(geno.exp), value=T)
    m.col <- grep("\\.m",names(geno.exp), value=T)
    m.counts <-  rowSums(geno.exp[,m.col,with=F])
    # select ASE input when total ase counts are above threshold and 
    A <- which(m.counts>=ase)
    if(nrow(geno.exp[A,][abs(rsnp)==1,])<n){
        return("Not enough individuals with ASE counts")
    } else {
        DT <- geno.exp[A,]
        DT[,n:=rowSums(DT[,n.col,with=F])]
        DT[, m:=m.counts[A]]
        DT[, p:=n/m][rsnp==-1, p:=1-p]
        DT[,g:=abs(rsnp)]
        ## calculate z-statistic
        if(sum(DT$g==1) >= (nrow(DT)-2)){ ## all or almost all hets, p.null=0.5
            DT <- DT[g==1,]
            zeta=(mean(DT$p)-0.5)/(sqrt(var(DT$p)/nrow(DT)))
        } else {
            DT[, het:=ifelse(g==1,"Yes","No")]
            N=DT[,.N,by=het]
            mean.p=DT[, mean(p), by=het]
            var.p=DT[,var(p), by=het]
            zeta <- diff(mean.p$V1)/sqrt(sum(var.p$V1/N$N))
        }
        cor.yg <- cor(abs(geno.exp$rsnp), geno.exp$y)
        return(v=c(cor.yg=cor.yg, z.p=zeta))
    }
}



#' function for getting p.ase from stan input
#'
#' This function allows you to compare the proportion of ase counts corrected by haplotype from the stan input for a snp-gene pair between hets and not hets and return the z.statistic
#' @param l list with stan input, output from stan.neg.beta.prob.eff
#' @keywords ase proportion stan input
#' @export
#' @return vector with z-statistic
#'
#' z.ase()

z.ase <- function(l){
    if(length(names(l))==4){## unname list, only the names for internal objects
        tmp <- sapply(1:length(l$n), function(x) (l$n[[x]]/l$gm[x,m]) %*% l$p[[x]])
        tmp[which(l$gm$g.ase==-1)] <-  1 - tmp[which(l$gm$g.ase==-1)] ## correct by hap
        DT <- data.table(g=abs(l$gm[,g.ase]), p= tmp)
    } else {
        tmp <- sapply(1:length(l[[1]]$n), function(x) (l[[1]]$n[[x]]/l[[1]]$gm[x,m]) %*% l[[1]]$p[[x]])
        tmp[which(l[[1]]$gm$g.ase==-1)] <-  1 - tmp[which(l[[1]]$gm$g.ase==-1)] ## correct by hap
        DT <- data.table(g=abs(l[[1]]$gm[,g.ase]), p= tmp)
    }
    ## calculate z-statistic
    if(sum(DT$g==1) >= (nrow(DT)-2)){ ## all or almost all hets, p.null=0.5
        DT <- DT[g==1,]
        zeta=(mean(DT$p)-0.5)/(sqrt(var(DT$p)/nrow(DT)))
    } else {
        DT[, het:=ifelse(g==1,"Yes","No")]
        N=DT[,.N,by=het]
        mean.p=DT[, mean(p), by=het]
        var.p=DT[,var(p), by=het]
        zeta <- diff(mean.p$V1)/sqrt(sum(var.p$V1/N$N))
    }
    
    return(zeta)
}

#' function for getting p.ase from stan input fixed, input is a list with 14 elements
#'
#' This function allows you to compare the proportion of ase counts corrected by haplotype from the stan input for a snp-gene pair between hets and not hets and return the z.statistic when the input is fixed haplotype
#' @param l list with stan input, output from fixhap.eff
#' @keywords ase proportion stan input
#' @export
#' @return vector with z-statistic
#'
#' z.fix.ase()

z.fix.ase <- function(l){
    if(length(names(l))==14){## unname list, only the names for internal objects
        tmp <- l$n/l$m
        tmp[which(l$gase==-1)] <-  1 - tmp[which(l$gase==-1)] ## correct by hap
        DT <- data.table(g=abs(l$gase), p= tmp)
    } else {
        tmp <- l[[1]]$n/l[[1]]$m
        tmp[which(l[[1]]$gase==-1)] <-  1 - tmp[which(l[[1]]$gase==-1)] ## correct by hap
        DT <- data.table(g=abs(l[[1]]$gase), p= tmp)
    }
    ## calculate z-statistic
    if(sum(DT$g==1) >= (nrow(DT)-2)){ ## all or almost all hets, p.null=0.5
        DT <- DT[g==1,]
        zeta=(mean(DT$p)-0.5)/(sqrt(var(DT$p)/nrow(DT)))
    } else {
        DT[, het:=ifelse(g==1,"Yes","No")]
        N=DT[,.N,by=het]
        mean.p=DT[, mean(p), by=het]
        var.p=DT[,var(p), by=het]
        zeta <- diff(mean.p$V1)/sqrt(sum(var.p$V1/N$N))
    }
    
    return(zeta)
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
    for(i in DT1$gene_id){
        log.c <- log(unname(unlist(DT2[gene_id==i,-"gene_id"])))
        g <- unname(unlist(DT1[gene_id==i,grep("_GT", names(DT1),value=T), with=F]))
    
        fit <- lm(log.c ~ g)

        ret[[i]] <- data.table(summary(fit)$coeff)[-1,][,SNP:=DT1[gene_id==i,id]][,gene_id:=i]
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


#' get SNPs with 95% CI for a param does not contain the null, to filter tag snps to further analyse
#'
#' This function allows you to select stan summaries with 95%CI not containing the null for a parameter from a list of stan summaries
#' @param a list of stan summaries, named with snp or no name
#' @param y name of parameter to extract
#' @param z value for the null hypothesis
#' @keywords non-null credible intervals
#' @export
#' @return DT with mean effect size, min, max, CI and  whether the null is within the CI
#' stan.no.null()

stan.no.null <- function(a,y="bj", z=0){
    cis.x <- stan.to.plot(a,y=y)
    if(!is.null(names(a))){
        cis.x[, SNP:=names(a)]
        }
    ## add whether the CI contains the null hypothesis, 
    cis.x[,bj.null.CI:=0][z>get(names(cis.x)[2]) & z<get(names(cis.x)[3]), bj.null.CI:=1]
    return(cis.x)
}

#' convert "X|Y" format into a matrix X Y as in hap format
#'
#' This function allows you to convert from vcf haplotypes to hap format
#' @param DT data table with vcf to convert, first 5 cols chrom, pos, id, ref, alt, GT with sample_GT
#' @keywords convert vcf to hap
#' @export
#' @return matrix with cols individuals (2 cols per ind) row snp, genotype per chromosome 0 or 1
#' vcf2hap()

vcf2hap <- function(DT){
    ind <- grep("_GT$",names(DT))
    DT[, id:= paste(POS,REF,ALT, sep=":")]
    tmp.l <- lapply(ind,function(i) {
        tmp <- Reduce(rbind,strsplit(DT[[i]],split="|"))
        if(!is.matrix(tmp)){
            tmp <- tmp[-2]
            tmp <- matrix(as.numeric(tmp), nrow=1)
        }  else {
            tmp <- tmp[,-2]
            tmp <- apply(tmp,2,as.numeric)
        }
        })
        mat <- Reduce(cbind,tmp.l)
    rownames(mat) <- DT$id
    return(mat)
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


#' Get correlation of haps from reference panel, to input in tag from GUESSFM, from Chris
#' 
#' calculates correlation by cols 
#' @param x matrix (example haplotypes for the snps, cols snps, rows samples)
#' @keywords correlation haplotypes reference panel
#' @export
#' @return matrix of correlatons
#' cor2 
cor2 <- function (x) {
  SD.x=apply(x,2, sd)
  if(any(SD.x==0)) stop("For some snps SD=0, remove and re-run")
  tmp = 1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE))
  return(tmp)
}
    

##' Derive tag SNPs for a SnpMatrix object using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters,
##' then cuts the tree at 1-tag.threshold
##' @param X matrix of haplotypes, each rows is an observed hap, cols snps
##' @param tag.threshold threshold to cut tree, default=0.99
##' @param quiet if FALSE (default), show progress messages
##' @param method method used for heirarchical clustering.  See hclust for options.
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @export
##' tag.noGT

tag.noGT <- function(X,tag.threshold=0.99, quiet=FALSE,method="single") {
  
  r2 <- (cor2(X))^2
  tmp=1-r2
  D <- stats::as.dist(1-r2)
  hc <- stats::hclust(D, method=method)
  clusters <- cutree(hc, h=1-tag.threshold)
  
  snps.use <- names(clusters)[!duplicated(clusters)]
  groups <- split(names(clusters),clusters)
  
  ## now process each group, picking best tag
  n <- sapply(groups,length)
  names(groups)[n==1] <- unlist(groups[n==1])
  for(i in which(n>1)) {
    g <- groups[[i]]
    ##     cat(i,g,"\n")
    ##     print(r2[g,g])
    a <- apply(r2[g,g],1,mean)
    names(groups)[i] <- g[ which.max(a) ]
  }
  groups <- new("groups",groups,tags=names(groups))
  
  ## check
  r2 <- r2[tags(groups),tags(groups)]
  diag(r2) <- 0
  ##   if(max(r2)==1) 
  ##     stop("max r2 still 1!")
  if(!quiet)
    message("max r2 is now",max(r2),"\n")
  return(as(groups,"tags"))
  
}

group.tags <- function(tags, keep) {
  groups <- tags[ names(tags) %in% keep ]
  groups <- split(names(groups), groups)
}


#####################################################################################

##' Calculate variance of expected genotype from a list with each element an individual and for each individual a vector of probabilities for GT=0,1 or 2.
##'
##' @title var.eg
##' @param l list, each element is for one individual, for each ind a vector with p(G)
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @export
##' var.eg

var.eg <- function(l){
    var_eg <-c()
    for(i in names(l)){
        tmp <- sapply(l[[i]]$NB$p.g, function(k) sum(as.numeric(names(k))*k))
        var_eg <- c(var_eg,var(tmp))
   
    }
    names(var_eg) <- names(l)
    return(var_eg)
}


        
##' Calculates variance of genotype for a group of snps based on reference panel
##'
##' @param file1 path to legend file
##' @param file2 path to hap file
##' @param x vector with snps coded as pos:ref:alt
##' @return named vector with var(G)
##' @export
##' var.rp

var.rp <- function(file1,file2,x){
    cw <- as.numeric(gsub(":.*","",x))
    cw <- c(min(cw),max(cw))
    tmp <- haps.range(file1,file2,cw)
    tmp <- tmp[x,]
    ## get GT from hap
    tmp2 <- sapply(seq(1,ncol(tmp),2), function(i) rowSums(tmp[,i:(i+1)]))
    ## get var per snp
    var.snp <- apply(tmp2,1,var)
    return(var.snp)
}


#' wrap functions to calculate var(E(G)), var(G), r2=var(E(G))/var(G) and abs(`log2(aFC)_mean.ngt`-`log2(aFC)_mean.gt`) into a data.frame
#'
#' This function allows you to to calculate var(E(G)), var(G), r2=var(E(G))/var(G) and abs(`log2(aFC)_mean.ngt`-`log2(aFC)_mean.gt`) into a data.frame
#' @param genes, character vector of genes to select inputs from, defaults to NULL to use all genes run by model
#' @param path charcater vector with path to files with stan.input
#' @param pattern character vector with pattern to match
#' @param noGT data table with Gene_id and tags to extract gene/snps pairs
#' @param le.file path to legend file for specific chromosome
#' @param hap.file path to haplotype file for specific chromosome
#' @keywords stan inputvar(E(G))
#' @export
#' @return input data table noGT with new cols: var(E(G)), var(G) and r2
#' var.e()

var.e <- function(genes=NULL,path,pattern,noGT,le.file,hap.file){
    ## stan input
    if(is.null(genes)){
        tmp <- lapply(list.files(path, pattern=pattern, full.names=T), readRDS)
        genes <- names(tmp) <- gsub(pattern, "",gsub(".*ENS", "ENS", list.files(path, pattern=pattern, full.names=T)))
    } else {
        tmp <- lapply(genes, function(i) readRDS(list.files(path=path, pattern=paste0(i,pattern), full.names=T)))
        names(tmp) <- genes
    }
    ## Extract snps per gene from noGT
    snps <- lapply(genes, function(i) which(names(tmp[[i]]) %in% noGT[Gene_id==i,tag.ngt]))
    names(snps) <- genes
    
    ngt.tgs <- mapply(function(x,y) x[y], tmp , snps, SIMPLIFY = FALSE)
    
    ## Calculate Var(E(G)) for each gene-snp pair
    var.e1 <-lapply(ngt.tgs, var.eg)                             
                  
    ## transform to DT
    var.e1 <- rbindlist(lapply(var.e1, function(i) data.table(snp=names(i), var.exp.g=i)), idcol="Gene_id")

    ## Calculate var(G)
    var.g <- var.rp(file1=le.file, file2=hap.file, x=unique(var.e1$snp))

    ## make new cols in noGT

    DT <- merge(noGT, data.table(snp=names(var.g),var.rp=var.g), by.x="tag.ngt",by.y="snp", all.x=TRUE)

    DT <- merge(DT,var.e1, by.x=c("Gene_id", "tag.ngt"), by.y=c("Gene_id","snp"), all.x=T)

    DT[,r2:=var.exp.g/var.rp][,abs.dif.log2.aFCg.ngt:=abs(`log2(aFC)_mean.ngt`-`log2(aFC)_mean.gt`)]

    
    return(unique(DT))
}

#' wrap function to calculate r2=var(E(G))/var(G) to use for info input on Btrecase.rna.noGT.R
#'
#' This function allows you to to calculate select snps with r2=var(E(G))/var(G) above specific threshold
#' @param stan.noGT input list for one gene
#' @param rp.r matrix with reference panel haplotype info, cols individuals, rows snps
#' @param info, cut-off to remove snps, when r2<info, remove
#' @keywords stan inputvar(E(G))
#' @export
#' @return named vector with r2 and names snp_id, r2 above threshold
#' info.cut()

info.cut <- function(stan.noGT, rp.r, info){
    
    ## Calculate Var(E(G)) for each snp 
    var.e1 <-var.eg(stan.noGT)
                  
    ## transform to DT
    #var.e1 <- rbindlist(lapply(var.e1, function(i) data.table(snp=names(i), ar.exp.g=i)), idcol="Gene_id")

    ## Calculate var(G)

    ## get GT from reference panel haps
    tmp2 <- sapply(seq(1,ncol(rp.r),2), function(i) rowSums(rp.r[,i:(i+1), drop=FALSE]))
    if(!is.matrix(tmp2)) tmp2 <- matrix(tmp2,  nrow=1,dimnames=list(unique(names(tmp2)), NULL))
    ## get var per snp
    var.g <- apply(tmp2,1,var)
    
    ## r2
    r2 <- var.e1/var.g

    r2 <- r2[r2>=info]
    
    return(r2)
}


########################## Genotyping from RNA #######################

#' Compare RNA with DNA genotyping
#'
#' This functions returns a data table with number of variants that were called both in RNA and DNA
#' @param rna, data table with rna , cols POS, id, others permitted
#' @param dna data table with dna cols, at least POS, and id
#' @param variable character vector with name of the col in rna to stratify by, used cummulative
#' @param range if variable is numeric, vector with values to use as cut-offs 
#' @keywords genyting common variants DNA RNA
#' @export
#' @return data table
#' var.r.d()

var.r.d <- function(rna,dna,variable=NULL,range=NULL){
    if(!is.null(variable)) {
        if(is.numeric(rna[[variable]])){
            tmp <- rbindlist(lapply(range, function(i) data.table(variable=i,vars_rna=nrow(rna[get(variable)<=i,]) , common_DNA_pos= length(which(rna[get(variable)<=i,POS] %in% dna$POS)), common_DNA_id= length(which(rna[get(variable)<=i,id] %in% dna$id)) ) ))
            variable <- paste0("cum.",variable)
        } else {
            tmp <- rbindlist(lapply(unique(rna[[variable]]), function(i) data.table(variable=i,vars_rna=nrow(rna[get(variable)==i,]) , common_DNA_pos= length(which(rna[get(variable)==i,POS] %in% dna$POS)), common_DNA_id= length(which(rna[get(variable)==i,id] %in% dna$id)) ) ))

        }
        setnames(tmp,"variable",variable)
    } else {
        tmp <- data.table(vars_rna=nrow(rna) , common_DNA_pos= length(which(rna[,POS] %in% dna$POS)), common_DNA_id= length(which(rna[,id] %in% dna$id)) )
    }
    tmp[, common_id_perc:=common_DNA_id * 100 / vars_rna]
    return(tmp)
}

#' Select common samples from dna and rna and look for concordance, adapted from Cinc/various.R
#'
#' This function allows you to: select variants that are in the same position and same REF and ALT in DNA and RNA; Add column Concordance for genotype; Add column counting # of genotyping discrepancies between RNA and DNA. RNA unphased, DNA phased
#' @param dna data table with DNA info, one chromosome
#' @param rna data table with RNA data (all chrs)
#' @param chr chromosome to subset
#' @param prefix whether sample names are different in RNA and DNA, if different provide prefix for RNA (without "_GT")
#' @keywords common variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype.
#' convert ()
convert2 <- function(dna,rna, chr, prefix=NULL){
    DT <- merge(rna[CHROM==chr,][,CHROM:=as.numeric(CHROM)], dna, by=c("CHROM","POS","REF","ALT"), suffixes=c("_RNA","_DNA"))

    dna_GT<-grep("GT", names(dna), value=TRUE)
    if(!is.null(prefix)){
        rna_GT<- paste0(prefix,"_GT")
    } else {
        prefix <- dna_GT  ## sample
        rna_GT <- paste0(dna_GT,"_RNA")  ## GTs
        dna_GT <- paste0(dna_GT,"_DNA")
    }  
    
    
    #cols<-paste0(cols_pre, "_GT")
for(i in 1:length(prefix)){
	## DT[,paste0(prefix[i],"_Concordance"):="Discordant"]
	## DT[get(paste0(rna_GT)) == "./.", paste0(prefix[i],"_Concordance"):=NA]
	## DT[get(dna_GT[i])==get(rna_GT[i]), paste0(prefix[i],"_Concordance"):="Concordant"]	
    DT[,paste0(prefix[i],"_N_errors"):=0]##[get(paste0(prefix[i],"_Concordance"))=="Concordant", paste0(prefix[i],"_N_errors"):=0]
    DT[get(rna_GT[i]) == "./.", paste0(prefix[i],"_N_errors"):=NA]
    DT[ get(dna_GT[i]) == "0|0" & (get(rna_GT[i]) == "1/0" | get(rna_GT[i]) == "0/1"), paste0(prefix[i],"_N_errors"):=1]
    DT[ get(dna_GT[i]) == "0|0" & get(rna_GT[i]) == "1/1", paste0(prefix[i],"_N_errors"):=2]
    DT[ get(dna_GT[i]) == "1|1" & get(rna_GT[i]) == "0/0", paste0(prefix[i],"_N_errors"):=2]  
    DT[ get(dna_GT[i]) == "1|1"  & (get(rna_GT[i]) == "1/0" | get(rna_GT[i]) == "0/1"), paste0(prefix[i],"_N_errors"):=1]   
    DT[ (get(dna_GT[i]) == "0|1" | get(dna_GT[i]) == "1|0" ) & (get(rna_GT[i]) == "1/1" | get(rna_GT[i]) == "0/0"), paste0(prefix[i],"_N_errors"):=1]
     
}
    return(DT)

}    


 #' get N-errors per sample for a chr, based on Cinciantti.various.R
#'
#' This function allows you to: summarize errors across samples per chr
#' @param dna_rna_10 data table with DNA and RNA info, one chromosome, output from convert2 function
#' @param format if not NULL gives output formatted for table
#' @param col variable to stratify results on (default null).
#' @keywords error calling variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' N_error2 ()
N_error2 <- function(dna_rna_10,format=NULL,col=NULL) {
        sum_errors <- N_error2_sub(dna_rna_10)
        temp1 <- N_error2_sub2(sum_errors,format)
      if(is.null(col)==TRUE){  
        return(temp1)
    } else {
        u <- unique(unlist(dna_rna_10[,col, with=F]))
        temp <- lapply(u, function(i) dna_rna_10[get(col)==i,])
        names(temp) <- u
        sum_errors_l <- lapply(temp,N_error2_sub)
        temp2 <- lapply(sum_errors_l, N_error2_sub2,format)
        ## add temp1 for summary
        temp2[["all"]] <- temp1
        temp2 <- rbindlist(temp2,idcol=col)
        setkeyv(temp2,col)
        return(temp2)
    }
}

#' subfunction 1 for get N-errors per sample for a chr, based on Cincinatti/various
#'
#' This function allows you to: summarize errors across samples per chr
#' @param dna_rna_10 data table with DNA and RNA info, one chromosome
#' @keywords common variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' N_error2_sub ()
N_error2_sub <- function(dna_rna_10){
    cols <- grep("GT_RNA",names(dna_rna_10), value=T) ## RNA samples
    DT<-data.table(Errors=c(0,1,2,NA))
    samples <- gsub("_RNA","",cols) ## prefix for _N_errors
 
    for (i in seq_along(samples)) {
        
        if(length(unique(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))]))==4){
            DT[,samples[i]:=as.numeric(data.table(table(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))], useNA="always"))$N)]
        } else{
            w <- which(DT$Errors %in% unique(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))]))
            DT[,samples[i]:=0]
            DT[w,samples[i]:=as.numeric(data.table(table(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))], useNA="always"))$N)]
        }
    }
    return(DT)

}
        
#' subfunction 2 for get N-errors per sample for a chr, modified from Cincinatti/various.R
#'
#' This function allows you to: summarize errors across samples per chr
#' @param DT data table with RNA errors relative to DNA
#' @param format returns table with formatted Mean (SD)
#' @keywords common variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' N_error2_sub2 ()

N_error2_sub2 <- function(DT, format=NULL) {
sum_errors<-data.table(Errors=c(0,1,2,NA), Mean=round(apply(DT[,2:ncol(DT)],1,mean)), SD=round(apply(DT[,2:ncol(DT)],1,sd)))
sum_errors[1:3,Pct:=round(Mean*100/sum(Mean[1:3]),2)]


    if(is.null(format)==TRUE){
        return(sum_errors)
        } else {
    sum_errors[,`Mean (SD)`:=paste0(Mean," (",SD,")")]
    sum_errors[,Mean:=NULL][,SD:=NULL]
    setcolorder(sum_errors, c("Errors","Mean (SD)","Pct"))
    return(sum_errors)
    
        }
}    


#' Get the proportion of concordance of RNA to DNA of particular GT by DP per sample
#'
#' This function allows you get the proportion of concordance of RNA to DNA of particular GT by DP per sample
#' @param DT data table with RNA and DNA genotypes, RNA depths
#' @param GT character vector with GT of DNA to consider ("homo", "het")
#' @param x value to group DP to visualise better concordance
#' @keywords concordance variants RNA DNA
#' @export
#' @return data table containing proportion of concordant GT in RNA by DP
#' err.dp ()

err.dp <- function(DT, GT=c("het","homo"), x=150) {
    ## samples prefix
    gt.dna <- grep("GT_DNA",names(DT), value=T)
    samp <- gsub("_GT_DNA","", gt.dna)
    gt <- c("1|0","0|1") ## dna het
    gt.rna <-c("1/0","0/1")  ## error in rna
    if(GT=="homo"){
        gt <- c("0|0","1|1")
        gt.rna <- c("0/0","1/1")
    }
    
    gt.e <- lapply(seq_along(samp), function(i) {
        ## concordance b DP
        dt <- DT[(get(gt.dna[i])==gt[1] | get(gt.dna[i])==gt[2]) &
                 (get(paste0(samp[i],"_GT_RNA"))==gt.rna[1] | get(paste0(samp[i],"_GT_RNA"))==gt.rna[2]), .N, by=get(paste0(samp[i],"_DP"))]

        ## total hets excluding missing values 
        dt2 <- DT[(get(gt.dna[i])==gt[1] | get(gt.dna[i])==gt[2]) &
                  get(paste0(samp[i],"_GT_RNA"))!="./.", .N, by=get(paste0(samp[i],"_DP"))]
        
        dt3 <- merge(dt,dt2,by="get",all.y=T,sort=T,suffixes=c(paste0(".ok.",i),paste0(".all.",i)))
        dt3[is.na(dt3)] <- 0
        return(dt3)
    })
    ## merge gt.e datatables to get concordance across samples   
    all <- Reduce(function(a,b) merge(a, b, by = "get", all = TRUE), gt.e)
    ## get results by DP
    conc <- rowSums(all[,grep("N.ok", names(all)),with=F], na.rm=T)
    all2 <- rowSums(all[,grep("N.all", names(all)),with=F], na.rm=T)
    ## group DP>=150
    p=data.table(p=c(conc[which(all$get %in% 1:x)], sum(conc[which(all$get %in% (x+1):nrow(all))]))/c(all2[which(all$get %in% 1:x)], sum(all2[which(all$get %in% (x+1):nrow(all))])), DP=all$get[1:(x+1)])
   return(p)
}

#' Get the DP across samples when GT is called
#'
#' This function allows you get the DP for each sample all combined for GT that were called
#' @param DT data table with RNA and DNA genotypes, RNA depths
#' @keywords concordance variants RNA DNA
#' @export
#' @return vector with DP for each sample and each GT called
#' all.dp ()

all.dp <- function(DT){
    ## samples prefix
    gt.rna <- grep("GT_RNA",names(DT), value=T)
    gt.dna <- grep("GT_DNA",names(DT), value=T)
    samp <- gsub("_GT_DNA","", gt.dna)
    dp <- paste0(samp,"_DP")
    tmp <- rbindlist(lapply(seq_along(samp), function(i) DT[get(gt.rna[i])!="./.",dp[i],with=F]))
    tmp2 <- Reduce("as.vector",tmp)
    return(tmp2)
}


#' Select GT for snps based on RNA by sample DP
#'
#' This function converts the RNA GT to ./. if DP <= value, for each sample
#' @param DT data table with RNA and DNA genotypes, RNA depths
#' @param DP cut-off for sample DP to consider GT, defaults=10
#' @keywords variant selection RNA 
#' @export
#' @return data table with corrected GT based on sample DP
#' gt_dp ()

gt_dp <- function(DT,x=10){
    gt.rna <- grep("GT_RNA",names(DT), value=T)
    dp <- grep("_DP", names(DT), value=T)
    tmp <- copy(DT)
    ## replace GT
    for( i in seq_along(gt.rna)){
        tmp[get(dp[i]) < x, gt.rna[i]:="./."]
    }
    ## remove N_errors col as they are not correct after recoding
    tmp[,grep("N_errors",names(tmp),value=T):=NULL]
    return(tmp)
}



#' Test of proportions for fsnps
#'
#' This function tests equality of the proportion of hets for fsnps in 2 populations, sample and reference panel
#' @param DT data table with fsnps and GT for samples (f.ase), cols for GT must end with "_GT"
#' @param matrix with rownames fsnps and col haplotypes of fsnps
#' @param gene character with gene_id
#' @keywords test hets proportion
#' @export
#' @return matrix with rownames fsnps and colnames OR (odds ratio) and pvalue for Fisher's exact test
#' prop_het ()

prop_het <- function(f.ase, rp.f,gene){
    ## get sample GT, count hets and totals
    sam.GT <- f.ase[,grep("_GT$",names(f.ase),value=T)]
    sam.rec <- rec_mytrecase_rSNPs(f.ase$POS,f.ase)
    sam.GT <- sam.rec[,grep("_GT$",names(f.ase),value=T), with=F]
    sam.hets <- apply(sam.GT,1, function(i) sum(abs(i)==1, na.rm=T))
    sam.all <- apply(sam.GT,1, function(i) sum(!is.na(i)))
    names(sam.hets) <- names(sam.all) <- f.ase$id
    sam.nohets <- sam.all-sam.hets

    ## get GT from reference panel haps, count hets and total=nrow(tmp)
    tmp <- sapply(seq(1,ncol(rp.f),2), function(i) rowSums(rp.f[,i:(i+1), drop=F]))
    if(!is.matrix(tmp)) tmp <- matrix(tmp, nrow=1,dimnames=list(unique(names(tmp)), NULL))
    
    ref.hets <- apply(tmp,1, function(i) sum(abs(i)==1))

    ## run fisher's exact test for each fsnp, takes hets and no-hets to work out odds ratio
    mat <- rbind(sam.hets[names(ref.hets)], sam.nohets[names(ref.hets)], ref.hets, ncol(tmp)-ref.hets)
    fish <- apply(mat,2, function(i){
        f=fisher.test(matrix(i,nrow=2),alternative="two.sided")
        v=c(f$estimate, f$p.value)
        names(v) <- c("OR", "pvalue")
        return(v)
    })

    fish <- data.table(t(fish), keep.rownames=T)
    fish[,gene_id:=gene]
    setnames(fish, "rn", "fsnp")
    return(fish)
    
}

