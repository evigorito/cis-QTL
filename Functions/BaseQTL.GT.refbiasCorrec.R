suppressMessages(library(ggplot2))
suppressMessages(library(MASS))
library(parallel)
suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


source('/home/ev250/Bayesian_inf/trecase/Functions/aux_btrecase_GT.R')


#' Run BaseQTL with known rsnp GT, optional refbias correction
#'
#' This function allows you to run BaseQTL for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
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
#' baseqtl.gt.refbias()

baseqtl.gt.refbias <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,u.esnps=NULL, gene.coord,vcf,le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9, out=".", prefix=NULL, model=c("both","trecase","trec"), stan.model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.refBias2.stan' ,
                                stan.trec='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.eff2.stan',
                                prob=NULL, prior=NULL, ex.fsnp=NULL, AI_estimate=NULL, pretotalReads=100) {
    
    ## prepare inputs
    base.in <- btrecase.gt.refbias.in(gene=gene,
                                      chr=chr,
                                      snps=snps,
                                      counts.f=counts.f,
                                      covariates=covariates,
                                      e.snps,
                                      u.esnps,
                                      gene.coord,
                                      vcf,
                                      le.file,
                                      h.file,
                                      population,
                                      nhets,
                                      min.ase,
                                      min.ase.het,
                                      tag.threshold,
                                      out,
                                      prefix,
                                      model,
                                      stan.model,
                                      stan.trec,
                                      prob,
                                      prior,
                                      ex.fsnp,
                                      AI_estimate,
                                      pretotalReads)
    
    if(is.character(base.in)) stop(base.in)

        ## get inputs
        if(any(names(base.in) == "trecase")) { ## proceed with full model
            stan.in2 <- base.in$trecase$trecase
            ASE.hets <- base.in$trecase$ASE.hets
            eaf.t <-base.in$trecase$eaf.t
            probs <- base.in$trecase$probs
            nhets <- base.in$trecase$nhets
            nfsnps <- base.in$trecase$nfsnps
            r.tag <- base.in$trecase$r.tag

            
            print(paste("Running full model", stan.model))

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
            names(stan.full) <- names(stan.in2)
            full.sum <- stan.bt(x=stan.full,y= NULL, rtag=r.tag, model="trec-ase", nhets=nhets, ASE.het=ASE.hets,gene=gene, EAF=eaf.t, nfsnps=nfsnps, probs=probs )

            ## add min AI_post
            if(!is.null(AI_estimate)) full.sum[, min_AI:= base.in$trecase$minAI]

        }

        if(any(names(base.in) == "neg")){
            in.neg <- base.in$neg$neg
            eaf.t <- base.in$neg$eaf.t
            r.tag <- base.in$neg$r.tag
            probs <- base.in$neg$probs
            nhets <- base.in$neg$nhets

            if(file.exists(stan.trec)){
                mod2 <- stan_model(stan.trec)
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
                names(stan.neg) <- names(in.neg)

                neg.sum <- stan.bt(x=stan.neg, y=NULL,rtag=r.tag,model="trec", nhets=nhets,gene=gene, EAF=eaf.t, nfsnps="NA", probs=probs)

            }
        }

        if( (exists("full.sum") & exists("neg.sum")) | exists("neg.sum") ) {
            if(exists("full.sum")){
                
                neg.sum <- rbind(full.sum, neg.sum, fill=TRUE)
            }
            
            if(!is.null(prefix)){             
                
                write.table(neg.sum, paste0(out,"/",prefix,".stan.summary.txt"), row.names=FALSE)
               
            } else {
               
                write.table(neg.sum,paste0(out,"/",gene,".stan.summary.txt"), row.names=FALSE)
               
            }
        }
        


        if(exists("full.sum") & !exists("neg.sum")) {
            
            if(!is.null(prefix)){             
                
                write.table(full.sum, paste0(out,"/",prefix,".stan.summary.txt"), row.names=FALSE)
                
            } else {
                
                write.table(full.sum,paste0(out,"/",gene,".stan.summary.txt"), row.names=FALSE)
                
            }
        }

}
            
            
             
        
                
                
                
        
