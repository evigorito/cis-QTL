
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

source('/home/ev250/Bayesian_inf/trecase/Functions/aux_btrecase.R')

##library(devtools)

##install_github("chr1swallace/GUESSFM", ref="groups")
library(GUESSFM)


#' Run Btrecase with known rsnp GT, optional refbias correction, preparing inputs more efficiently
#'
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the chromosome where the gene is
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to get EAF for rsnp: AFR AMR EAS EUR SAS ALL, defaults to EUR
#' @param nhets minimun number of het individuals in order to run the minumn model (NB only), defaults to NULL
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.het minimun number of het individuals with the minimun of ASE counts in order to run the ASE side of the model, defaults to NULL
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to NULL
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving tables, if NULL gene_id.eqtl will be used
#' @param model whether to run trecase, trec or both, defaults to both
#' @param stan.model full name for stan model to run when using ASE side (with or w/o ref bias correction)
#' @param stan.trec full name for stan model to run trec only
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.95
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Deffaults to NULL
#' @param ex.fsnp, if character: vector with pos:ref:alt for fsnps to exclude,  defaults to NULL
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, defaults to NULL
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @keywords bayesian trecase known genotype regulatory snp
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from trecase model. Saves rds file with the names of the fSNPS used for inference (prefix/GT.fsnps.with.counts.rds). Saves rds file with the inputs to run the model for QC purposes (prefix/GT.stan1.input.rds).
#' btrecase.gt.refbias()

btrecase.gt.refbias <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,u.esnps=NULL, gene.coord,vcf,le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9, out=".", prefix=NULL, model=c("both","trecase","trec"), stan.model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.refBias2.stan' ,
                                stan.trec='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.eff2.stan',
                                prob=NULL, prior=NULL, ex.fsnp=NULL, AI_estimate=NULL, pretotalReads=100) {
    
    ## check inputs and extract inputs for gene
    model <- match.arg(model)
    if(!is.null(prior)){
        if(class(prior) != "list") stop("prior argument must be a list")
        if(any(!names(prior) %in% c("mean", "sd", "mix"))) stop("prior argument must be a named 'mean' and 'sd'")
        if(length(unique(sapply(prior, length))) !=1) stop("mean and sd for prior argument must have the same length")
    }
    ingene <- aux.in1(gene,
                      chr,
                      snps,
                      counts.f,
                      covariates,
                      e.snps,
                      u.esnps,
                      gene.coord,
                      vcf,
                      sample.file=NULL, ## only for noGT
                      le.file,
                      h.file,
                      population=population,
                      maf=NULL,
                      nhets=nhets,
                      min.ase=min.ase,
                      min.ase.het=min.ase.het,
                      min.ase.snp=NULL,
                      min.ase.n=NULL,
                      tag.threshold=tag.threshold,
                      out=out,
                      model=model,
                      stan.model=stan.model,
                      prob=prob,
                      ex.fsnp=ex.fsnp,
                      AI_estimate=AI_estimate,
                      pretotalReads=pretotalReads)

    counts.g <- ingene$counts
    covariates <- ingene$covariates
    probs <- ingene$probs
    
    ## get rsnps and extract GT and ASE, remove non-informative snps (missing or hom in all samples)
    ## create dt to collect rsnps excluded from analysis

        
    rsnps.ex <- data.table(id=character(), reason=character())
    gcoord <- fread(gene.coord)
    
    if(is.numeric(snps)){     
        if("percentage_gc_content" %in% names(gcoord)){## newer version of gcoord from psoriasis snakefile rule geno_info using start, end, chrom, longest transcript length and GC percentage
            cis_window <- tryCatch({gcoord[gene_id==gene & chrom==chr,.(start,end)] + c(-snps, snps)},
                                   error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        } else { ## old version fp gcoord
            cis_window <- tryCatch({cl_coord(file=gene.coord,chr,gene=gene,cw=snps)}, error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        }
        if(is.character(cis_window)) stop(cis_window)
        if(is.list(cis_window)) cis_window <- unlist(cis_window)
        gt.as <- vcf_w(vcf,chr, st=cis_window["start"], end=cis_window["end"], exclude="yes")
        if(is.character(gt.as)) stop(print(gt.as))
        rsnps.ex <- gt.as$excluded
        gt.as <- gt.as$keep
        ## make sure gt.as doesn't have missing values or unphased GT   
        gt.qc <- vcf.gt.qc(gt.as)
        if(sum(gt.qc[c(2,4)])!=0) stop(cat("Invalid GT field for some snps or samples \n",paste(names(gt.qc), collapse=","), "\n",gt.qc, "\n","to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments"))
        rs <- copy(gt.as)
        
    } else {
        pos <- as.numeric(sapply(strsplit(snps, split=":"), `[[`,1))
        w <- which(!is.na(pos))
        if(!length(w)) stop(cat("Invalid format for snps ", snps[w]))
        ## get gene start and end, ciswindow=0
        
        if("percentage_gc_content" %in% names(gcoord)){## newer version of gcoord from psoriasis snakefile rule geno_info using start, end, chrom, longest transcript length and GC percentage
            st_end <- tryCatch({gcoord[gene_id==gene & chrom==chr,.(start,end)] + rep(0, 2)},
                               error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        } else {
            st_end=tryCatch({cl_coord(gene.coord,chr,gene=gene,cw=0)}, error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        }

        if(is.character(st_end)) stop(st_end)
        if(is.list(st_end)) st_end <- unlist(st_end)
        ## construct cis-window with snps, making sure to include the whole gene
        m <- min(pos) < st_end[1]
        M <- max(pos) > st_end[2]
        cis_window <- ifelse(m, min(pos), st_end[1])
        cis_window <- c(cis_window, ifelse(M, max(pos), st_end[2]))
        gt.as <- vcf_w(vcf,chr, cis_window["start"], cis_window["end"], exclude = "yes")
        if(is.character(gt.as)) stop(print("snps not found in vcf"))
        rsnps.ex <- gt.as$excluded[id %in% snps,]
        gt.as <- gt.as$keep
        rs <- gt.as[id %in% snps,]
        if(!nrow(rs)) stop("Missing GT or homozygous snps in all samples")        
        ## make sure rsnps doesn't have missing values or unphased GT   
        gt.qc <- vcf.gt.qc(rs)
        if(sum(gt.qc[c(2,4)])!=0) stop(cat("Invalid GT field for some snps or samples \n",paste(names(gt.qc), collapse=","), "\n",gt.qc, "\n","to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments"))
    }
    ## further process of rsnps
    rs <- rs[,grep("_AS",names(rs),value=T):=NULL]       
    ## recodehel to 0,1,-1,2 scale
    rec.rs <- rec_mytrecase_rSNPs(x=rs$POS, y=rs)
    if(tag.threshold!="no") {
        ## Group rsnps by r2, recode rec.rs for input in tags function from GUESSFM
        if(nrow(rs)==1) stop("Only one regulatory snp to test, please set tag.threshold='no' \n Cannot cluster one snp only")
        re.guess <- rec.guess(DT=rec.rs)
        x <- as(re.guess-1, "SnpMatrix")
        rtag <- GUESSFM::tag(X=x,tag.threshold=tag.threshold)
        ## save rtag as data.table 
        dt <- data.table(Gene_id=gene,tag=tags(rtag), SNP=rtag@.Data)
        if(!is.null(prefix)){
            write.table(dt,paste0(out,"/",prefix,".tags.lookup.txt"), row.names=FALSE)
        } else {
            write.table(dt,paste0(out,"/",gene,".eqtl.tags.lookup.txt"), row.names=FALSE)
        }
        
        ## restrict rsnp to tag snps
        rec.rs <- rec.rs[id %in% unique(tags(rtag)),]
    }
    r.tag <- switch(is.numeric(tag.threshold),"yes") ## to output results after stan, when tag.threshold is char, returns NULL
    ## help with hets
    GT.aux <- rec.rs[,grep("_GT",names(rec.rs)),with=F] ## to make easier calculation of correlations, etc
    
    ## counts number of hets per rsnp
    rec.rs[, nhet:=apply(GT.aux ,1, function(i) sum(abs(i)==1))] 
    ## remove snps with less than min hets
    w <- rec.rs[nhet<nhets, which = TRUE]
    rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs[w,id], reason=rep(paste("rsnp with less than", nhets ,"het ind."), length(w))))
    rec.rs <- rec.rs[!w,]
    if(nrow(rec.rs)==0) stop(cat("No rsnp with at least", nhets ,"het ind."))

    
    if(model=="trecase" | model=="both") {
        
        ## get fSNPs (feature snps or exonic snps)
        fsnps <- tryCatch({fread(cmd=paste("grep", gene, e.snps))}, error=function(e) {paste("No entry for gene",gene, "in",e.snps)})
        
        if(isTRUE(grep("No entry",fsnps)==1) | nrow(fsnps)==0 ) { ## no fsnps
            cat(fsnps, "\n No fSNPS, analysis will be done with total gene counts only")
            rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs$id, reason=paste("No entry for gene",gene, "in",e.snps)))
        } else {     
            fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
            fsnps <- fsnps[id %in% gt.as$id]
            ## exclude fsnps if required
            if(!is.null(ex.fsnp)){
                fsnps <- fsnps[!id %in% ex.fsnp,]
            }

            ## Only use fSNPs with AI_estimates for SNPs with sufficient reads
            
            if(!is.null(AI_estimate)) {
                ai <- fread(AI_estimate)    
                ai <- ai[CHROM == chr & Total_pre >=pretotalReads & Keep == "yes" & AI_post != 0,]
                ai[, id:=paste(POS,REF,ALT, sep = ":")]

                fsnps <- fsnps[ id %in% ai$id, ]

                if(nrow(fsnps) == 0) {
                    cat("No fSNPs with allelic imbalance estimates \n analysis will be done with total gene counts only")
                    rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs$id, reason="No fSNPs with AI estimates"))
                }
            } else {

                ai <- NULL
            }
            
            ##
            
            if(is.character(snps)) { ## analysis on pre-specified set of rsnps
                gt.as <- gt.as[id %in% c(rec.rs$id,fsnps$id),] ## snps and fsnps only
                ## make sure gt.as (here for fsnps) doesn't have missing values or unphased GT   
                gt.qc <- vcf.gt.qc(gt.as)
                if(sum(gt.qc[c(2,4)])!=0) stop(cat("Invalid GT field for some exonic snps and samples \n",paste(names(gt.qc), collapse=","), "\n",gt.qc, "\n","to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments \n", "if all exonic snps are removed the analysis will be done with total counts only"))           
            }
            
            ## get info from reference panel 
            ## matrix for reference panel haps
            rp <- haps.range(file1=le.file,file2=h.file,cis_window,maf=0)
            ## keep record of rsnps not in reference panel
            w <- which(!rec.rs$id %in% rownames(rp))
            rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs$id[w],reason=rep("No rsnp in reference panel", length(w))))
            rec.rs2 <- rec.rs[!w,] ## for full model, only rsnps in ref panel
            w <- which(rownames(rp)%in%  rec.rs2$id)
            if(length(w)!=0) {  ## proceed with full model, neg only for rsnps not in ref panel           
                w=which(rownames(rp)%in% gt.as$id)
                rp <- rp[w,, drop=FALSE]
                ## make sure to select fsnps from reference panel
                f.ase <- gt.as[id %in% fsnps$id,]
                f.ase <- f.ase[id %in% rownames(rp),]
                
                if(nrow(f.ase)==0) {
                    print("No fsnps ref panel")
                    rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs2$id, reason="No fsnps in reference panel"))
                } else {
                    
                    ## get counts per hap for fsnp
                    counts.ase <- tot.ase_counts(x=f.ase)

                    if(!is.null(u.esnps)){
                        counts.ase <- aux.in2(gene, u.esnps, counts.ase, ex.fsnp, ai)
                        
                    }
                    
                    if(!is.character(counts.ase)) { ## enough fsnps for ASE

                        ## check if sufficient indiviuals with ASE counts to proceed

                        if(!is.character(filt.fsnp(counts.ase,min.ase.snp=0))) { ## Enough ind with ase counts
                            
                            print(paste("Effective number of fSNPs:", nrow(f.ase)))                        
                            
                            ## select rsnps in ref panel for full model
                            rs.full <- rs[id %in% rec.rs2$id, which(names(rs) %in% names(rec.rs2)), with=F]
                            
###################  run stan full model #######################
###### prepare stan inputs
                            ## reference panel fsnps and rsnps:
                            rp.f <- rp[f.ase$id,,drop=FALSE]
                            rp.r <- rp[rs.full$id,,drop=FALSE]

                            ## get haplotype counts from fsnps

                            ## add AI_estimate per sample to stan inputs
                            if(!is.null(ai)) ai <- ai[id %in% f.ase$id,.(id, NREF_post, NALT_post,Total_post, AI_post)]

                            print("Preparing stan inputs")

                            ## order fsnps in couts.ase as in f.ase 

                            counts.ase=counts.ase[,unlist(lapply(f.ase$id, function(i) grep(i, colnames(counts.ase), value = TRUE)))]

                            stan.f <- stan.fsnp.noGT.eff(rp.f,f.ase,counts.ase, NB="no", min.ase, min.ase.het, ai=ai)

                            
                            if(!is.character(stan.f)) { ## proceed with full model
                                
                                counts <- unlist(counts.g) ## to avoid repeating in mclapply
                                
                                stan.in1<-mclapply(rs.full$id, function(i) stan.trecase.eff2(counts, rp.1r=rp.r[i,,drop=FALSE], rp.f, f.ase, rs.hap=rs.full[id==i,], rec.rsnp=rec.rs2[id==i,], stan.f,min.ase, min.ase.het))
                                
                                names(stan.in1) <- rs.full$id

                                ## remove rsnps with insufficient ase, min.ase.hets (haplotypes with fsnps not compatible with reference panel)
                                w <- sapply(stan.in1, is.character)
                                if(any(w)){
                                    rsnps.ex <- rbind(rsnps.ex,data.table(id=rs.full$id[w], reason="Not enough het individuals with sufficient ase counts or Hap not in reference panel"))   
                                    
                                    stan.in1 <- stan.in1[!w]
                                }
                                if(length(stan.in1) >= 1) {
                                                                       
                                    ## save stan.in1 to qc
                                    if(is.null(prefix)) {
                                        saveRDS(stan.in1, paste0(out,"/",gene,".GT.stan1.input.rds"))
                                        saveRDS(colnames(counts.ase)[colSums(counts.ase)>5], paste0(out,"/",gene,".GT.fsnps.with.counts.rds"))
                                    } else{
                                        saveRDS(stan.in1, paste0(out,"/",prefix,".GT.stan1.input.rds"))
                                        saveRDS(colnames(counts.ase)[colSums(counts.ase)>5], paste0(out,"/",prefix,".GT.fsnps.with.counts.rds"))
                                    }
                                    
                                    print(paste("Running full model", stan.model))
                                    
                                    if(is.null(prior)){
                                        stan.in2 <- lapply(stan.in1, function(i) in.neg.beta.prob.eff2(i, covar=covariates))
                                    } else {
                                        stan.in2 <- lapply(stan.in1, function(i) {
                                            l <- in.neg.beta.prob.eff2(i, covar=covariates)
                                            ## add prior
                                            ## if k=1 I need to make as.array, otherwise stan takes vector as real
                                            k=unique(sapply(prior,length))
                                            l[['k']] =unique(sapply(prior,length))
                                            if(k==1) {
                                                l[['aveP']]=as.array(prior$mean)
                                                l[['sdP']]=as.array(prior$sd)
                                                ## log of mixing proportions to avoid calculation in stan
                                                l[['mixP']]=as.array(log(prior$mix))
                                            } else {
                                                l[['aveP']]=prior$mean
                                                l[['sdP']]=prior$sd
                                                ## log of mixing proportions to avoid calculation in stan
                                                l[['mixP']]=log(prior$mix)
                                            }
                                            return(l)
                                        }      
                                        )
                                    }

                                    mod <- stan_model(stan.model)
                                  
                                    stan.full <- mclapply(stan.in2, function(i) {
                                        post <- sampling(mod,data=i, cores=1, refresh=0 , pars='bj')
                                        s <- summary(post, pars='bj', use_cache=F, probs=probs)$summary
                                        ##samp <- sampling(mod,data=i, cores=1, refresh=0, pars="bj")
                                        ## s <- summary(samp, pars='bj', use_cache=F, probs=probs)$summary
                                        e <- rstan::extract(post, pars="bj")
                                        ## calculate proportion of post <0
                                        post.neg <- sum(e$bj<0)/length(e$bj)
                                        ## add to s
                                        s <- cbind(s['bj',, drop=F], post.prop.neg=post.neg)
                                        unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                                        return(s)
                                    })

                                    names(stan.full) <- names(stan.in1)
                                    ## count number of hets with sufficient ase counts to input in final report
                                    ASE.hets <- sapply(stan.in1, function(i) nrow(i$gm[abs(g.ase)==1,]))
                                    ## get eaf for tags run in model
                                    eaf.t <- snp.eaf(le.file,names(stan.full),population)
                                    full.sum <- stan.bt(x=stan.full,y= NULL, rtag=r.tag, model="trec-ase", nhets=rec.rs2[id%in% names(stan.full),nhet], ASE.het=ASE.hets,gene=gene, EAF=eaf.t, nfsnps=nchar(unlist(strsplit(colnames(stan.f$n[[1]]), ","))[1]), probs=probs )

                                    ## add min AI_post
                                    if(!is.null(AI_estimate)) full.sum[, min_AI:= min(ai$AI_post)]
                                    
                                } ## closing from failed stan.in1
                            } else {  ## all rsnps will fail because stan.f failed
                                rsnps.ex <- rbind(rsnps.ex,data.table(id=rs.full$id,  reason=stan.f))
                                
                            }## closing from failed stan.f
                        } else { ## not Enough inds with ase counts
                       
                            w <- which(!rec.rs2$id  %in% rsnps.ex$id)
                            if(length(w)) rsnps.ex <- rbind(rsnps.ex, data.table(id= rec.rs2$id[w], reason ="No enough inds with ASE counts"))
                        }
                        
                    } else {  ##  no unique fSNPs, enough fsnps for ASE

                        rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs2$id, reason="Not unique fsnps in gene"))
                        
                    }  ## closing no unique fSNPs
                   
                } ## closing from no fsnps in ref panel
               
            } ## closing from rsnps in ref panel
            
        } ## closing from fsnps in gene
        
                
        if(model == "trecase") {
            if(nrow(rsnps.ex)){
                if(!is.null(prefix)) {             
                    write.table(rsnps.ex,paste0(out,"/",prefix,".trecase.excluded.rsnps.txt"), row.names=FALSE)
                    
                } else {
                    write.table(rsnps.ex,paste0(out,"/",gene,".trecase.excluded.rsnps.txt"), row.names=FALSE)
                }
            }
            if(exists("full.sum")){
                
                if(!is.null(prefix)) {             
                    
                    write.table(full.sum, paste0(out,"/",prefix,".stan.summary.txt"), row.names=FALSE)
                    ##saveRDS(list(full=stan.post), paste0(out,"/",prefix,".stan.post.RDS"))
                    
                    
                } else {
                    
                    write.table(full.sum,paste0(out,"/",gene,".stan.summary.txt"), row.names=FALSE)
                    ##saveRDS(list(full=stan.post), paste0(out,"/",gene,".stan.post.RDS"))
                }
                return(full.sum)  
            } else {
                return("No rsnps were run with trecase")
            }
        }

        
    } else { ## when model==trec I add all rec.rs snps to prepare input

        rsnps.ex <- rbind(rsnps.ex,data.table(id=rec.rs$id, reason="requested model trec"))

    }
    
    ## run neg.binom only 
    ## prepare input and run
    if(nrow(rsnps.ex)>0) {
        
        ex <- c("Missing or homo GT all samples", "rsnp with less than","below cut-off for qtest")
        w <- unlist(sapply(ex, function(i) grep(i,rsnps.ex$reason)))
        id.keep <- unique(rsnps.ex[!w,id])
        
        if(length(id.keep)) {
            if(is.null(prior)) {
                in.neg <- mclapply(id.keep, function(i) input.neg.only.bj(DT1=counts.g, DT2=rec.rs[id ==i,] ,covar=covariates))
            } else {
                in.neg <- mclapply(id.keep, function(i) {
                    l <- input.neg.only.bj(DT1=counts.g, DT2=rec.rs[id ==i,] ,covar=covariates)
                    ## add prior
                    ## if k=1 I need to make as.array, otherwise stan takes vector as real
                    l[['k']] =unique(sapply(prior,length))
                    k=unique(sapply(prior,length))
                    if(k==1) {
                        l[['aveP']]=as.array(prior$mean)
                        l[['sdP']]=as.array(prior$sd)
                        ## log of mixing proportions to avoid calculation in stan
                        l[['mixP']]=as.array(log(prior$mix))
                    } else {
                        l[['aveP']]=prior$mean
                        l[['sdP']]=prior$sd
                        ## log of mixing proportions to avoid calculation in stan
                        l[['mixP']]=log(prior$mix)
                    }
                    
                    
                    return(l)
                })
            }
            
           
            mod2 <- stan_model(stan.trec)
            names(in.neg) <- id.keep
           
            print(paste("Running NB model", stan.trec))
            ## get full posterior
            stan.neg <-  mclapply(in.neg, function (i) {
                samp <- sampling(mod2,data=i, cores=1, refresh=0, pars="bj")
                s <- summary(samp, pars='bj', use_cache=F, probs=probs)$summary
                ## extract params from posterior
                e <- rstan::extract(samp, pars="bj")
                ## calculate proportion of post <0
                post.neg <- sum(e$bj<0)/length(e$bj)
                ## add to s
                s <- cbind(s['bj',, drop=F], post.prop.neg=post.neg)
                
                unload.ddl(mod2) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                return(s)
            })
            names(stan.neg) <- id.keep
        
            ## get eaf for tags run in model
            eaf.t <- snp.eaf(le.file,names(stan.neg),population)
            ## if tag not in ref panel (possible with neg binom side)
            if(length(id.keep) > nrow(eaf.t)){
                eaf.t <- rbind(eaf.t, data.table(snp=id.keep[which(!id.keep %in%  eaf.t$snp)], eaf=NA))
                ## sort as in id.keep
                eaf.t <- eaf.t[order(match(snp, id.keep))]
            }

            
            neg.sum <- stan.bt(x=stan.neg, y=NULL,rtag=r.tag,model="trec", nhets=rec.rs[id %in% id.keep,nhet],gene=gene, EAF=eaf.t, nfsnps="NA", probs=probs)
            if( exists("full.sum")){
                neg.sum <- rbind(full.sum, neg.sum, fill=TRUE)
                ## setorder(neg.sum,`log2_aFC_null`,-`log2_aFC_d.aux`, -model)
               
                
            }
            if(!is.null(prefix)){             
                write.table(rsnps.ex,paste0(out,"/",prefix,".trecase.excluded.rsnps.txt"), row.names=FALSE)
                write.table(neg.sum, paste0(out,"/",prefix,".stan.summary.txt"), row.names=FALSE)
               
            } else {
                write.table(rsnps.ex,paste0(out,"/",gene,".eqtl.trecase.excluded.rsnps.txt"), row.names=FALSE)
                write.table(neg.sum,paste0(out,"/",gene,".stan.summary.txt"), row.names=FALSE)
               
            }          
            return(neg.sum)
        }

    } else {
        if(model == "trec"){
            
            return("No rsnps run using trec")
            
        }
        if(exists("full.sum")){
            if(!is.null(prefix)) {             
                write.table(full.sum, paste0(out,"/",prefix,".stan.summary.txt"), row.names=FALSE)
                write.table(rsnps.ex,paste0(out,"/",prefix,".trecase.excluded.rsnps.txt"), row.names=FALSE)
            } else {
                write.table(full.sum,paste0(out,"/",gene,".stan.summary.txt"), row.names=FALSE)
                write.table(rsnps.ex,paste0(out,"/",gene,".eqtl.trecase.excluded.rsnps.txt"), row.names=FALSE)
                
            }
            return(full.sum)  
        } else {         
            return("No snps run using trecase or trec")
        }    
    }
    if(!exists("full.sum")){
        if(nrow(rsnps.ex)){
         if(!is.null(prefix)){             
                write.table(rsnps.ex,paste0(out,"/",prefix,".trecase.excluded.rsnps.txt"), row.names=FALSE)                           
            } else {
                write.table(rsnps.ex,paste0(out,"/",gene,".eqtl.trecase.excluded.rsnps.txt"), row.names=FALSE)
                              
            }
         }
         return("No snps were run")
    }
}





