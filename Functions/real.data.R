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
    ##snp.st <- paste0("/",snps)
    ##sed.mult <- paste(snp.st, collapse="/= ;") ## quicker but I loose track of snps if some arent on ref panel
    
    line <- lapply(1:length(snps), function(i) paste0("zcat ", file1, " | sed -n  '/", snps[i], "/='"))  #sed a little bit quicker than awk and grep
    
    nlines <- lapply(line, function(i) as.numeric(system(i, intern=TRUE))) ## slow but need to kepp track of snps in case some are not in reference panel 
    names(nlines) <- snps
    
    ## sort before reading
    nlines <- sort(unlist(nlines))

    if(length(nlines)==0){
        return("no snps in reference panel")

    } else {
    
    
        ## get haplotypes for fsnps in ref panel, extract from hap.gz file
        ## first line in legend file is headings, need to substract 1 to match hap.gz file
        
        p.nlines <- paste0(nlines-1, "p")
        s.nlines <- nlines-1

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
                    if(sum(vl[(w[i-1]+1):w[i]])==2){ ## the new stretch of 1's in w[i] starts one snp after the last stretch of ones
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
        rownames(mat) <- names(nlines)
        return(mat)
    }
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


#' function for filtering input for stan
#'
#' This function allows you to test wether a rsnp has enough het inds with a certain number of ASE counts
#' @param geno.exp data table for an element of output list tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5
#' @keywords filtering stan input
#' @export
#' @return data table with sample name, GT rsnp, total counts, n and m counts for each fSNP, same data table outputted as an element of tot.ase_counts
#'
#' filt.rsnp()

filt.rsnp <- function(geno.exp,ase=5, n=5){
   
    n.col <- grep("\\.n",names(geno.exp), value=T)
    m.col <- grep("\\.m",names(geno.exp), value=T)
    m.counts <-  rowSums(geno.exp[,m.col,with=F])
    # select ASE input when total ase counts are above threshold and 
    A <- which(m.counts>=ase)
    if(nrow(geno.exp[A,][abs(rsnp)==1,])<n){
        return("Not enough individuals with ASE counts")
    } else {
        return(geno.exp)
    }
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


    

    
    
    
    
