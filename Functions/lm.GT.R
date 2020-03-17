
suppressMessages(library(ggplot2))
suppressMessages(library(MASS))
library(parallel)
suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

##library(devtools)

##install_github("chr1swallace/GUESSFM", ref="groups")
library(GUESSFM)

#' Run lm with known rsnp GT 
#'
#' This function allows you to run lm for one gene and multiple pre-selected snps.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the chromosome where the gene is
#' @param gt file with exrated genotypes from vcf, when using this option set vcf to null, default is to use vcf
#' @param nhets minimun number of het individuals in order to run the model, defaults to 5
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving tables, if NULL gene_id.eqtl will be used
#' @param pval whether to return pvalue, defaults to no return (ci only)
#' @keywords lm known genotype regulatory snp
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from trecase model.
#' lm.gt()

lm.gt <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,gene.coord, vcf, gt=NULL,nhets=5,tag.threshold=.9, out=".", prefix=NULL, pval=NULL) {
    
    ## check inputs
    files <- c(counts.f,gene.coord)
    if(!is.null(vcf)) {
        files <- c(files, vcf)
    }else{
        files <- c(files, gt)
    }
    
    w <- !file.exists(files)
    if(any(w)) stop(paste("invalid file names ", paste(files[w], collapse= ", ")))
    if(!(tag.threshold >=0 & tag.threshold <=1 | tag.threshold =="no")) stop("invalid  'tag.threshold' argument")
    if(!(dir.exists(out))) stop("invalid 'out' argument")
    num <- nhets
    names(num) <-  "nhets"
    w <- !is.numeric(num)
    if(any(w)) stop(paste("invalid arguments:",names(num)))
    if(!(covariates==1 | file.exists(covariates))) stop("invalid path for covariates file")
    
    ## Extract inputs for gene

    ## get counts 
    counts.g <- fread(cmd=paste("grep -e gene_id -e ",gene,counts.f), header=TRUE)
    counts.g <- counts.g[,2:ncol(counts.g),with=F] ## removes gene_id
    if(nrow(counts.g)==0) stop("Gene id is not found in count matrix")
    if(covariates !=1) {
        covariates <- readRDS(covariates)
        if(!(is.matrix(covariates))) stop("Covariates file is not a matrix")
        if(gene %in% rownames(covariates)){
            ## covariates by gene
            covariates <- covariates[gene,, drop=F]
            covariates <- t(covariates)
        } else {
            if(ncol(covariates) >= ncol(counts.g)){
                stop("Either the number of covariates is >= number of samples or the gene is not in covariates")
            }
        }      
        
        if(!(nrow(covariates)==ncol(counts.g)))  stop(cat("Number of individuals in covariates matrix:", nrow(covariates), "\n Number of individuals in gene counts", ncol(counts.g), "\n please adjust"))
        ## scale
        covariates <- apply(covariates,2,scale,center=TRUE,scale=TRUE)
        rownames(covariates)=names(counts.g)
    }
    
    ## get rsnps and extract GT and ASE, remove non-informative snps (missing or hom in all samples)
    ## create dt to collect rsnps excluded from analysis
    rsnps.ex <- data.table(id=character(), reason=character())
    gcoord <- fread(gene.coord)
    if(is.numeric(snps)) {
         if("percentage_gc_content" %in% names(gcoord)){## newer version of gcoord from psoriasis snakefile rule geno_info using start, end, chrom, longest transcript length and GC percentage
            cis_window <- tryCatch({gcoord[gene_id==gene & chrom==chr,.(start,end)] + c(-snps, snps)},
                                   error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        } else { ## old version fp gcoord
            cis_window <- tryCatch({cl_coord(file=gene.coord,chr,gene=gene,cw=snps)}, error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        }
        
         if(is.character(cis_window)) stop(cis_window)
         if(is.null(gt)){
             gt.as <- vcf_w(vcf,chr, st=cis_window["start"], end=cis_window["end"], exclude="yes")
             if(is.character(gt.as)) stop(print(gt.as))
             rsnps.ex <- gt.as$excluded
             gt.as <- gt.as$keep
             ## make sure gt.as doesn't have missing values or unphased GT   
             gt.qc <- vcf.gt.qc(gt.as)
             if(sum(gt.qc[c(2,4)])!=0) stop(cat("Invalid GT field for some snps or samples \n",paste(names(gt.qc), collapse=","), "\n",gt.qc, "\n","to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments"))
             rs <- copy(gt.as)
         } else {
             gt <- readRDS(gt)
             gt[, id:=paste(POS,REF,ALT, sep=":")]
             rec.rs <- gt[POS>=cis_window$start & POS <= cis_window$end,]
         }
         
        
    } else {
        pos <- as.numeric(sapply(strsplit(snps, split=":"), `[[`,1))
        w <- which(!is.na(pos))
        if(!length(w)) stop(cat("Invalid format for snps ", snps[w]))
        if(!is.null(vcf)){
            ## get gene start and end, ciswindow=0
            if("percentage_gc_content" %in% names(gcoord)){## newer version of gcoord from psoriasis snakefile rule geno_info using start, end, chrom, longest transcript length and GC percentage
                st_end <- tryCatch({gcoord[gene_id==gene & chrom==chr,.(start,end)] + rep(0, 2)},
                                   error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
            } else {
                st_end=tryCatch({cl_coord(gene.coord,chr,gene=gene,cw=0)}, error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
            }

            ## construct cis-window with snps, making sure to include the whole gene
            #m <- min(pos) < st_end[1]
            #M <- max(pos) > st_end[2]
            cis_window <- setNames(c(min(pos, unlist(st_end)), max(pos,unlist(st_end))), c("start", "end"))
            gt.as <- vcf_w(vcf,chr, cis_window["start"], cis_window["end"], exclude = "yes")
            if(is.character(gt.as)) stop(print("snps not found in vcf"))
            rsnps.ex <- gt.as$excluded[id %in% snps,]
            gt.as <- gt.as$keep
            rs <- gt.as[id %in% snps,]
            if(!nrow(rs)) stop("Missing GT or homozygous snps in all samples")        
            ## make sure rsnps doesn't have missing values or unphased GT   
            gt.qc <- vcf.gt.qc(rs)
            if(sum(gt.qc[c(2,4)])!=0) stop(cat("Invalid GT field for some snps or samples \n",paste(names(gt.qc), collapse=","), "\n",gt.qc, "\n","to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments"))
        } else {
            gt <- readRDS(gt)
            gt[, id:=paste(POS,REF,ALT, sep=":")]
            rec.rs <- gt[id %in% snps,]
        }
    }
    ## further process of rsnps
    if(exists("rs")){
        rs <- rs[,grep("_AS",names(rs),value=T):=NULL]       
    ## recode to 0,1,-1,2 scale
        rec.rs <- rec_mytrecase_rSNPs(x=rs$POS, y=rs)
    }
    
    if(tag.threshold!="no") {
        ## Group rsnps by r2, recode rec.rs for input in tags function from GUESSFM
        if(nrow(rs)==1) stop("Only one regulatory snp to test, please set tag.threshold='no' \n Cannot cluster one snp only")
        re.guess <- rec.guess(DT=rec.rs)
        x <- as(re.guess-1, "SnpMatrix")
        rtag <- GUESSFM::tag(X=x,tag.threshold=tag.threshold)
        ## save rtag as data.table 
        dt <- data.table(Gene_id=gene, tag=tags(rtag), SNP=rtag@.Data)
        if(!is.null(prefix)){
            write.table(dt,paste0(out,"/",prefix,".lm.tags.lookup.txt"), row.names=FALSE)
        } else {
            write.table(dt,paste0(out,"/",gene,".lm.eqtl.tags.lookup.txt"), row.names=FALSE)
        }
        
        ## restrict rsnp to tag snps
        rec.rs <- rec.rs[id %in% unique(tags(rtag)),]
    }
    r.tag <- switch(is.numeric(tag.threshold),"yes") ## to output results after lm, when tag.threshold is char, returns NULL
  
    ## help with hets and GT
    GT.aux <- rec.rs[,grep("_GT",names(rec.rs)),with=F] ## to make easier calculation of correlations, etc

    ## remove rsnps with same GT in all samples
    rem <- apply(GT.aux,1, function(i) length(unique(abs(i))) == 1)
    if(any(rem)){
        
        rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs[rem,id], reason=rep(paste("Same GT for rsnp in all samples."), length(rem))))
        rec.rs <- rec.rs[!rem,]
        GT.aux <- GT.aux[!rem,]
        
    }
    if(all(rem)) stop("Each rsnp have the same genotype in every sample")

    ## counts number of hets per rsnp
    rec.rs[, nhet:=apply(GT.aux ,1, function(i) sum(abs(i)==1))]

    ## remove snps with less than min hets
    w <- rec.rs[nhet<nhets, which = TRUE]
    if(length(w)) {
        rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs[w,id], reason=rep(paste("rsnp with less than", nhets ,"het ind."), length(w))))
        rec.rs <- rec.rs[!w,]
        if(nrow(rec.rs)==0) stop(cat("No rsnp with at least", nhets ,"het ind."))
    }
    if(nrow(rsnps.ex)){
        rsnps.ex[,Gene_id:=gene]
        setcolorder(rsnps.ex,c("Gene_id", "id","reason"))
            if(!is.null(prefix)) {             
                write.table(rsnps.ex,paste0(out,"/",prefix,".lm.excluded.rsnps.txt"), row.names=FALSE)
                
            } else {
                write.table(rsnps.ex,paste0(out,"/",gene,".lm.excluded.rsnps.txt"), row.names=FALSE)
            }
    }       
    ## run lm
    ## if any sample has 0 counts, add 1 count to all samples to run model
    if(any(counts.g==0)){
        counts <- unlist(counts.g + 1)
    } else {
        counts <- unlist(counts.g)
    }
    
    log.c <- log(counts)
    if(!is.null(pval)){
        lm.fit <- data.table(aFC_mean=numeric(), pvalue=numeric(), aFC_se=numeric(),  aFC_2.5=numeric(), aFC_97.5=numeric())
    } else {
        lm.fit <- data.table(aFC_mean=numeric(),aFC_se=numeric(), aFC_2.5=numeric(), aFC_97.5=numeric())
    }
    
    for(i in 1:nrow(rec.rs)){
        g <- abs(unlist(rec.rs[i, grep("_GT", names(rec.rs),value=T), with=F]))
        ## create formula
        ## get cols in covariates
        if(is.matrix(covariates)){
            col.cov <- sapply(1:ncol(covariates), function(i) paste0("covariates[,",i,"]"))
            fmla <- paste("log.c ~ g + ", paste(col.cov, collapse="+"))
        } else {
            fmla <- "log.c ~ g"
        }
             
        fit <- lm(as.formula(fmla))
        ## express beta genotype  as beta allelic imbalance in log2 scale
        beta <- coefficients(fit)['g']
        bai <- 2*beta/log(2)  ## allelic fold change is beta*2, convert to log2 scale
        ## use p-value to get cIs in log2 scale
        p=summary(fit)$coeff["g","Pr(>|t|)"]
        #z=qnorm(p/2,lower.tail=FALSE)
        #ci <- c(bai-1.96*abs(bai)/z , bai+1.96*abs(bai)/z)
        ########### working here ##############
        se <- 2*summary(fit)$coeff["g","Std. Error"]/log(2) ## 2* to account for bai instead of beta
        ci <-  c(bai-1.96*se , bai+1.96*se)
        if(!is.null(pval)){
            lm.fit <- rbind(lm.fit, data.table(aFC_mean=bai, pvalue=p,  aFC_se= se, aFC_2.5=ci[1], aFC_97.5=ci[2]))
        } else {
            lm.fit <- rbind(lm.fit, data.table(aFC_mean=bai, aFC_se= se, aFC_2.5=ci[1], aFC_97.5=ci[2]))
        }
        
    }

    lm.fit[, aFC_null:= ifelse(0>=aFC_2.5 & 0<= aFC_97.5, "yes", "no")]
    lm.fit[aFC_null=="no" & aFC_2.5 >0 ,aFC_d:= aFC_2.5][aFC_null=="no" & aFC_2.5 <0 ,aFC_d:= -aFC_97.5][aFC_null=="yes",aFC_d:=abs(aFC_2.5 - aFC_97.5)]
    lm.fit[, aFC_d.aux:=aFC_d][aFC_null=="yes", aFC_d.aux:=-aFC_d] ## to sort length CI in ascending order
    names(lm.fit) <- paste0('log2.', names(lm.fit))   
    if(!is.null(r.tag)){
        lm.fit[,tag:=rec.rs$id]
    } else {
        lm.fit[,SNP:=rec.rs$id]
    }
    lm.fit[, Gene_id:=gene]
    lm.fit[,model:="lm-log(counts)"]
    setorder(lm.fit,log2.aFC_null,-log2.aFC_d.aux)
    if(!is.null(r.tag)){
        setcolorder(lm.fit, c("Gene_id", "tag", grep("log2", names(lm.fit), value=T), "model"))
    } else {
        setcolorder(lm.fit, c("Gene_id", "SNP", grep("log2", names(lm.fit), value=T), "model"))
    }
    
    
    if(!is.null(prefix)){             
        write.table(lm.fit, paste0(out,"/",prefix,".lm.summary.txt"), row.names=FALSE)
    } else {
        write.table(lm.fit, paste0(out,"/",gene,".lm.summary.txt"), row.names=FALSE)
    }          
    return(lm.fit)
}











