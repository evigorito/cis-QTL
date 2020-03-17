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


##################################################################################################################


#' Run Btrecase with unknown rsnp GT and missing values for GT fsnps with optional  reference panel bias correction
#'
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. 
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R. If using gc correction (each gene diffrent value), the matrix has rownames= genes and cols=samples plus extra columns if other covariates are added. If only using lib size or gene independent covariates, rows are samples and columns are covariates. If no covariates, covariates =1, default
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the exonic snps for the chromosome where the gene is
#' @param sample.file sample file for the reference panel (sample description), to be used if ex.fsnp test is required and populatios is not the whole reference panel
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL, and for testing fSNPs to exclude (ex.fsnp below) if applicabloe, defaults to EUR
#' @param maf cut-off for maf, defaults to 0.05
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.snp minum number of ASE counts for a single snp to be considered, for a particular individual
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving files, if NULL gene_id.eqtl will be used
#' @param model full name  to stan model, defaults to /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.priors.refbias.stan
#' @param model.negonly stan model with neg only side, deafults to NULL
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.95
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Deffaults to NULL
#' @param ex.fsnp, if character: vector with pos:ref:alt for fsnps to exclude, if numeric  p-value cut-off for fSNPs to exclude,  defaults to NULL
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @keywords bayesian trecase unknown genotype regulatory snp reference panel bias
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from trecase model.
#' btrecase.nogt.rna.refbias()

btrecase.nogt.rna.refbias <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,u.esnps=NULL, gene.coord,vcf,sample.file=NULL, le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), maf=0.05, min.ase=5,min.ase.snp=5,min.ase.n=5,tag.threshold=.9, info=0.3, out=".", prefix=NULL, model="/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.noGT.rsnp.prior02.refbias2.stan", model.negonly=NULL, prob=NULL, prior=NULL, ex.fsnp=NULL, AI_estimate=NULL, pretotalReads=100) {
  
    inputs <- btrecase.nogt.rna.refbias.In(gene=gene,
                                           chr=chr,
                                           snps=snps,
                                           counts.f= counts.f,
                                           covariates=covariates,
                                           e.snps = e.snps,
                                           u.esnps=u.esnps,
                                           gene.coord=gene.coord,
                                           vcf=vcf,
                                           sample.file=sample.file,
                                           le.file=le.file,
                                           h.file=h.file,
                                           population=population,
                                           maf=maf,      
                                           min.ase=min.ase,
                                           min.ase.snp=min.ase.snp,
                                           min.ase.n=min.ase.n,
                                           tag.threshold=tag.threshold,
                                           info=info,
                                           out=out,
                                           model=NULL,
                                           prefix=prefix,
                                           prob=prob,
                                           ex.fsnp=ex.fsnp,                                                                                           AI_estimate=AI_estimate,
                                           pretotalReads=pretotalReads)

    if(!is.null(prior)){
        if(class(prior) != "list") stop("prior argument must be a list")
        if(any(!names(prior) %in% c("mean", "sd", "mix"))) stop("prior argument must be a named 'mean' and 'sd'")
        if(length(unique(sapply(prior, length))) !=1) stop("mean and sd for prior argument must have the same length")
    }

    
    if(is.character(inputs$inp)) stop(inputs$inp)

    stan.noGT2 <- inputs$inp$stan.noGT2
    if(!is.null(ex.fsnp) & is.numeric(ex.fsnp)) {
        het.f <- inputs$inp$het.f
        het.fall <- inputs$inp$het.all
    } else {
        het.f=NULL
    }
    
    nfsnps <- inputs$inp$nfsnps
    info.ok <- inputs$inp$info.ok
    probs <- inputs$probs
    r.tag <- inputs$r.tag
    
    

    if(!is.null(prior)){
        stan.noGT2 <- lapply(stan.noGT2, function(l){ 
            ## add prior
            l[['k']] =unique(sapply(prior,length))
            l[['aveP']]=prior$mean
            l[['sdP']]=prior$sd
            ## log of mixing proportions to avoid calculation in stan
            l[['mixP']]=log(prior$mix)
            return(l)
        }
        )
    }

    if(inputs$model== "full"){
        model2run <- model
    } else {
        if(!is.null(model.negonly)){
            ## neg only model =="NB"
            if(file.exists(model.negonly)){
                model2run <- model.negonly
            } else {
                stop("Wrong name for model.negonly file")
            }
        }
    }

    if(!is.null(model2run)) {
        cat("Running stan with ", model2run)  
        
        mod <- stan_model(model2run)
        model <- inputs$model #"BNB-ASE"
        
        stan.full <-  mclapply(1:length(stan.noGT2),
                               function (i) {
                                   post <- sampling(mod,data=stan.noGT2[[i]], cores=1, refresh=0 , pars='bj')
                                   s <- summary(post, use_cache=F, probs=probs)$summary
                                   ## extract params from posterior
                                   e <- rstan::extract(post, pars="bj")
                                   ## calculate proportion of post <0
                                   post.neg <- sum(e$bj<0)/length(e$bj)
                                   ## add to s
                                   s <- cbind(s['bj',, drop=F], post.prop.neg=post.neg)
                                   unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                                   return(s)
                               })
        names(stan.full) <- names(stan.noGT2)
       
        ## get maf for tags run in model
        maf.t <- snp.eaf(le.file,names(stan.noGT2),population)
        if(model == "NB"){
            nfsnps=NA
        }
        ## stan.bt creates PEP from post.prop.neg
        full.sum <- stan.bt(x=stan.full,y= "bj",rtag=r.tag,model=model ,gene=gene,EAF=maf.t, info=info.ok, nfsnps=nfsnps, min.pval=het.f, probs=probs)
        ## add min AI_post
        if(!is.null(AI_estimate) & model == "full"){
            full.sum[, min_AI:= inputs$inp$min_AI]
        }
        
        if(!is.null(prefix)) {             
            
            write.table(full.sum, paste0(out,"/",prefix,".noGT.stan.summary.txt"), row.names=FALSE)
            if(exists("het.fall")){
                write.table(het.fall, paste0(out,"/",prefix,".fsnps.het.fisher.test.txt"), row.names=FALSE)
            }
            
        } else {
            
            write.table(full.sum,paste0(out,"/",gene,".noGT.stan.summary.txt"), row.names=FALSE)
            if(exists("het.fall")){
                write.table(het.fall, paste0(out,"/",gene,".fsnps.het.fisher.test.txt"), row.names=FALSE)
            }
            
        }
        
        return(full.sum)  
    } else {
        return("No model run")
    }

}



