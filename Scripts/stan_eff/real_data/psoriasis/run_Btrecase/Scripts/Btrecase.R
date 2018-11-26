## log <- file(snakemake@output[[1]], open="wt")
## sink(log)
## sink(log, type="message")

## sink()

suppressMessages(library(ggplot2))
suppressMessages(library(MASS))
library(parallel)
suppressMessages(library(rstan))
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

library(GUESSFM)

#' Run Btrecase with unknown rsnp GT and missing values for GT fsnps
#'
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. 
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the exonic snps for the chromosome where the gene is
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL, defaults to EUR
#' @param maf cut-off for maf, defaults to 0.05
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.snp minum number of ASE counts for a single snp to be considered, for a particular individual
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param q.test whether to select snps to run based on a pre-filtering test, character string indicating "yes", "no" (default)
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving files, if NULL gene_id.eqtl will be used
#' @param model path to stan model, defaults to /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.priors.eff2.stan
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.95
#' @param ex.fsnp, character vector with pos:ref:alt for fsnps to exclude, defaults to NULL
#' @keywords bayesian trecase unknown genotype regulatory snp
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from trecase model.
#' btrecase.nogt.rna()

btrecase.nogt.rna <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,gene.coord,vcf,le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), maf=0.05, min.ase=5,min.ase.snp=5,min.ase.n=5,tag.threshold=.9,q.test="no", info=0.3, out=".", prefix=NULL, model="/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.eff2.stan", prob=NULL, ex.fsnp=NULL) {
  
    ## check inputs
    ##model <- match.arg(model)
    population <- match.arg(population)
    files <- c(counts.f,e.snps,gene.coord,vcf,le.file,h.file)
    w <- !file.exists(files)
    if(any(w)) stop(paste("invalid file names ", paste(files[w], collapse= ", ")))
    if(!(tag.threshold >=0 & tag.threshold <=1 | tag.threshold =="no")) stop("invalid  'tag.threshold' argument")
    if(!is.null(prob)) {
        if (prob < 0 | prob >1 ) stop("invalid 'prob' argument")
    }
    
    probs <- if(is.null(prob)) NULL else c((1-prob)/2, 0.25, 0.5, 0.75, (1+prob)/2) 
    na.q.test <- pmatch(q.test, c("yes","no"))
    if(is.na(na.q.test)) stop("invalid 'q.method' argument")
    if(!(dir.exists(out))) stop("invalid 'out' argument")
    num <- list(maf,min.ase,min.ase.snp, min.ase.n,info)
    names(num) <-  c("maf", "min.ase","min.ase.snp", "min.ase.n","info")
    w <- !sapply(num, is.numeric)
    if(any(w)) stop(paste("invalid arguments:",paste(names(num)[w], collapse=", ")))
    if(covariates !=1){
        if(!file.exists(covariates)) stop("invalid path for covariates file")
    }
    
    if(is.character(snps)){
        pos <- as.numeric(sapply(strsplit(snps, split=":"), `[[`,1))
        w <- which(!is.na(pos))
        if(!length(w)) stop(cat("Invalid format for snps ", snps[w]))
    }   

    ## Extract inputs for gene
    
    ## get counts  and covariates
    counts.g <- fread(paste("grep -e gene_id -e ",gene,counts.f), header=TRUE)
    if(nrow(counts.g)==0) stop("Gene id is not found in count matrix")
    counts.g <- counts.g[,2:ncol(counts.g),with=F] ## removes gene_id
    ## get covariates and scale
    if(covariates !=1){
        covariates <- t(readRDS(covariates))
        if(!(is.matrix(covariates))) stop("Covariates file is not a matrix")
        if(nrow(covariates)!=ncol(counts.g))  stop(cat("Number of individuals in covariates matrix:", nrow(covariates), "\n Number of individuals in gene counts", ncol(counts.g)-1, "\n please adjust"))
        ## scale
        covariates <- apply(covariates,2,scale,center=TRUE,scale=TRUE)
        rownames(covariates)=names(counts.g)
    }
    
    ## create dt to collect snps excluded from analysis
    snps.ex <- data.table(id=character(), reason=character())
    
    ## get fSNPs (feature snps or exonic snps)
    fsnps <- fread(paste("grep", gene, e.snps))

    ## exclude fsnps if required
    if(!nrow(fsnps)) stop(cat("No entry for gene ", gene, " in e.genes"))
   
    fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
    if(!is.null(ex.fsnp)){
        fsnps <- fsnps[!id %in% ex.fsnp,]
    }
    if(!nrow(fsnps)) stop(cat("No entry for gene ", gene, " in e.genes"))
    
    ## get fsnps and extract GT and ASE, remove non-informative snps (missing or hom in all samples)
    fsnps.window <- c(min(fsnps$V2),max(fsnps$V2))
    gt.as <- vcf_w(vcf,chr, st=fsnps.window[1], end=fsnps.window[2], exclude="yes")
    if(is.character(gt.as)) stop(print(gt.as))
    snps.ex <- gt.as$excluded[id %in% fsnps$id,]
    gt.as <- gt.as$keep[id %in% fsnps$id]
    if(!nrow(gt.as)) stop("No exonic snps in vcf")
    
    ## get ASE counts per individual
    c.ase <- tot.ase_counts(x=gt.as)
    
    ## Check if there are sufficient hets with enough ASE counts
    c.ase = filt.fsnp(c.ase,ase=min.ase,min.ase.snp,n=min.ase.n, rem=NULL)
    if(is.character(c.ase)) stop(c.ase)
    
  
    ## Extract haps for fsnps and rsnp from reference panel ########
    if(is.numeric(snps)) {
        cis_window <- tryCatch({cl_coord(file=gene.coord,chr,gene,cw=snps)}, error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        if(is.character(cis_window)) stop(cis_window)
        rp <- haps.range(file1=le.file,file2=h.file,cw=cis_window,population, maf)
        if(is.character(rp)) stop(rp)
        if(!nrow(rp)) stop("No snps extracted from the reference panel")
        rp.r=rp
    } else {
        pos <- as.numeric(sapply(strsplit(snps, split=":"), `[[`,1))
        w <- which(!is.na(pos))
        if(!length(w)) stop(cat("Invalid format for snps ", snps[w]))
        ## get gene start and end, ciswindow=0
        st_end=tryCatch({cl_coord(gene.coord,chr,gene,cw=0)}, error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        if(is.character(st_end)) stop(st_end)
        ## construct cis-window with snps, making sure to include the whole gene
        m <- min(pos) < st_end[1]
        M <- max(pos) > st_end[2]
        cis_window <- ifelse(m, min(pos), st_end[1])
        cis_window <- c(cis_window, ifelse(M, max(pos), st_end[2]))
        rp <- haps.range(file1=le.file,file2=h.file,cis_window,population, maf)
        if(!nrow(rp)) stop("No snps extracted from the reference panel")
        w=which(!snps %in% rownames(rp))
        if(length(w)==length(snps)) stop(cat("None of the snps ", snps, "are in the reference panel"))   
        if(length(w)) {
            rp.r=rp[snps[-w],,drop=FALSE] ## only rsnps requested and in ref panel
            snps.ex=rbind(snps.ex, data.table(id=snps[w], reason=rep("snp not in reference panel",length(w))))
        } else {
            rp.r <- rp[snps,,drop=FALSE]
        }
        

    }
    w=which(!fsnps$id %in% rownames(rp))
    if(length(w)==length(fsnps$id)) stop("None of the fsnps are in the reference panel")
    snps.ex=rbind(snps.ex, data.table(id=fsnps$id[w], reason=rep("exonic snp is not in reference panel with maf filtering", length(w))))
    
    ## only use fsnps in ref panel (filtered by maf)
    f.ase <- gt.as[id %in% rownames(rp),]
    
    rp.f=rp[rownames(rp) %in% f.ase$id, ,drop=FALSE] ## exonic snps to work with
    c.ase=c.ase[,unlist(lapply(f.ase$id, function(i) grep(i, colnames(c.ase), value = TRUE)))]
    ## Check if there are sufficient hets with enough ASE counts if removing fsnps not in reference panel
    if(nrow(rp.f) < nrow(fsnps)) {
        c.ase = filt.fsnp(c.ase,ase=min.ase,min.ase.snp,n=min.ase.n, rem=NULL)
        if(is.character(c.ase)) stop(c.ase)
    }
    print(paste("Effective number of exonic SNPs:", nrow(f.ase)))
    
    if(tag.threshold!="no") {
      ## Group rsnps by r2, use tags function from GUESSFM
      if(nrow(rp.r)==1) stop("Only one regulatory snp to test, please set tag.threshold='no' \n Cannot cluster one snp only")
      ## exclude snps with SD=0
      x=t(rp.r)
      SD.x=apply(x,2,sd)
      w=SD.x==0
      if(all(w)) stop("All SNPs have zero standard deviation, please correct maf cut-off and rerun")
      if(any(w)) {
        snps.ex=rbind(snps.ex, data.table(id=colnames(x)[w], reason=rep("Snp with zero standard deviation", sum(w))))
        x=x[,!w, drop=FALSE]
      }
      rtag <- tag.noGT(X=x,tag.threshold=tag.threshold)
      ## save rtag as data.table 
      dt <- data.table(Gene_id=gene, tag=tags(rtag), SNP=rtag@.Data)
      if(!is.null(prefix)){
        write.table(dt,paste0(out,"/",prefix,".noGT.tags.lookup.txt"), row.names=FALSE)
      } else {
        write.table(dt,paste0(out,"/",gene,".noGT.eqtl.tags.lookup.txt"), row.names=FALSE)
      }
      
      ## restrict rsnp to tag snps
      rp.r <- rp.r[unique(tags(rtag)),]
    }
    r.tag <- switch(is.numeric(tag.threshold),"yes") ## to output results after stan, when tag.threshold is char, returns NULL
    ## help with hets
 

##############################  Stan inputs ###############################
                    
    ## fsnps, pre-compute what I need because applies to every snp
    
    stan.f <- fsnp.prep2(rp.f, f.ase, c.ase , min.ase, min.ase.n)

    if(is.character(stan.f)){
        return(stan.f)
    }
            

    stan.noGT<-mclapply(1:nrow(rp.r), function(i) stan.trecase.rna.noGT.eff2(counts.g, rp.1r=rp.r[i,,drop=FALSE], rp.f, stan.f))
    names(stan.noGT) <- rownames(rp.r)
    
    ## remove rsnps when genotypes of fsnps are not compatible with reference panel)
    w <- sapply(stan.noGT, is.character)
    if(any(w)){
        snps.ex <- rbind(snps.ex,data.table(id=rownames(rp.r)[w], reason=unlist(stan.noGT[w])))   
        stan.noGT <- stan.noGT[!w]
        rp.r <- rp.r[names(stan.noGT),]
    }
    ## restrict rsnps to those with info over cut-off
    info.ok <- info.cut(stan.noGT,rp.r,info)
    ## remove snps below cut-off
    w <- !names(stan.noGT) %in% names(info.ok)
    if(any(w)){
        snps.ex <- rbind(snps.ex,data.table(id=names(stan.noGT)[w], reason="Below info cut-off"))   
        stan.noGT <- stan.noGT[!w]        
    }   
    
    if(nrow(snps.ex)){
        ## add gene and re-order
        snps.ex[, Gene_id:=gene]
        setcolorder(snps.ex,c("Gene_id", "id","reason"))
        if(!is.null(prefix)) {             
            write.table(snps.ex,paste0(out,"/",prefix,".noGT.trecase.excluded.snps.txt"), row.names=FALSE)
            
        } else {
            write.table(snps.ex,paste0(out,"/",gene,".noGT.trecase.excluded.snps.txt"), row.names=FALSE)
        }
        
    }
    ## prepare stan inputs
    if(length(stan.noGT) >= 1) {
        ## save inputs for qc
        ##saveRDS(stan.noGT, paste0(out,"/", gene, ".noGT.stan.in.rds"))
        stan.noGT2 <- lapply(stan.noGT, function(i) in.neg.beta.noGT.eff2(i, covar=covariates[names(i$NB$counts),, drop=F]))

        ## trec model if no ase, same input
        w <- sapply(stan.noGT, function(i) is.character(i$ase))
        if(any(w)){
            stan.trec <- stan.noGT2[w]
            stan.noGT2 <- stan.noGT2[!w]
        }
        
        if(length(stan.noGT2)){     
            mod <- stan_model(model)
            stan.full <-  mclapply(1:length(stan.noGT2),
                                   function (i) {
                                       if(is.null(prob)){
                                           s <- summary(sampling(mod,data=stan.noGT2[[i]], cores=1, refresh=0), pars='bj', use_cache=F)$summary
                                       } else {
                                           s <- summary(sampling(mod,data=stan.noGT2[[i]], cores=1, refresh=0), pars='bj', use_cache=F, probs=probs)$summary
                                       }
                                       
                                       unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                                       return(s)
                                   })
            names(stan.full) <- names(stan.noGT2)

            ## get maf for tags run in model
            maf.t <- snp.eaf(le.file,names(stan.noGT2),population)
            ## get test for het proportion for fsnps
            het.f <- suppressWarnings(prop_het(f.ase,rp.f,gene)) ## missing values converted to NA
            full.sum <- stan.bt(x=stan.full,y= "bj",rtag=r.tag,model="trec-ase",gene=gene,EAF=maf.t, info=info.ok, nfsnps=length(unique(Reduce(c, stan.f$f.comb))), min.pval=het.f, probs=probs)
           
        } 
        if(exists("stan.trec")){ ## trec model
 
            mod.neg.nogt <-  stan_model('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.noGT.rsnp.priors.eff2.stan')
            neg.only <- mclapply(1:length(stan.noGT2),
                                 function (i) {
                                     s <- summary(sampling(mod,data=stan.noGT2[[i]], cores=1, refresh=0))$summary
                                     unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                                     return(s)
                                 })
            names(neg.only) <- names(stan.trec)
            ## get maf for tags run in model
            maf.t <- snp.eaf(le.file,names(stan.trec),population)
            neg.sum <- stan.bt(x=neg.only,y= "bj",rtag=r.tag,model="trec",gene=gene,EAF=maf.t, info=info.ok, nfsnps="NA" )
           
        }
        
        if(exists("full.sum") & exists("neg.only")){
            full.sum <- rbind(full.sum, neg.sum)
            setorder(neg.sum,log2_aFC_null,-log2_aFC_d.aux, -model)
        } else {
            if(exists("neg.only")){
                full.sum=neg.only
            }
        } 
        
        if(!is.null(prefix)) {             
            
            write.table(full.sum, paste0(out,"/",prefix,".noGT.stan.summary.txt"), row.names=FALSE)
            if(exists("het.f")){
                write.table(het.f, paste0(out,"/",prefix,".fsnps.het.fisher.test.txt"), row.names=FALSE)
            }
            
        } else {
            
            write.table(full.sum,paste0(out,"/",gene,".noGT.stan.summary.txt"), row.names=FALSE)
            if(exists("het.f")){
                write.table(het.f, paste0(out,"/",gene,".fsnps.het.fisher.test.txt"), row.names=FALSE)
            }
            
        }
        
        return(full.sum)  
        
    } else {
        
        return("None of the snps  met the conditions to be run by trec-ase")
    }
    
}   


btrecase.nogt.rna(gene=snakemake@params[['gene']],
                  chr=as.numeric(snakemake@params[['chrom']]),
                  snps=as.numeric(snakemake@params[['snps']]),
                  counts.f=snakemake@input[['counts']],
                  #covariates=snakemake@input[['libsize']],
                  e.snps=snakemake@input[['eSNPs']],
                  gene.coord=snakemake@input[['genecoord']],
                  vcf=snakemake@input[['vcf']],
                  le.file=snakemake@input[['leRef']],
                  h.file=snakemake@input[['hapRef']],
                  population=snakemake@params[['pop']],
                  maf=as.numeric(snakemake@params[['maf']]),
                  min.ase=as.numeric(snakemake@params[['minAse']]),
                  min.ase.snp=as.numeric(snakemake@params[['minAseSnp']]),
                  min.ase.n=as.numeric(snakemake@params[['minAseN']]),
                  tag.threshold=as.numeric(snakemake@params[['tag']]),
                  q.test="no",
                  info=as.numeric(snakemake@params[['info']]),
                  out=snakemake@params[['out']],
                  prefix=NULL,
                  model="/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.eff2.stan",
                  prob=NULL,
                  ex.fsnp=NULL)


#e.snps='/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/fSNP/chr1.fSNP.genes.txt'
#gene="ENSG00000233645"
gene="ENSG00000000457"
counts.f <-  "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/Counts/Psoriasis_skin.txt"
vcf <- "/mrc-bsu/scratch/ev250/psoriasis/Btrecase/inputs/GT/chr1.ASE.Psoriasis_skin.vcf.gz"
chr <- 1
le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz'

h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz'

model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.eff2.stan'

population="EUR"
maf=0.05
nhets=5
min.ase=5
min.ase.snp=5
min.ase.n=5
tag.threshold=.9
q.test="no"
info=0.3
snps=5*10^5
