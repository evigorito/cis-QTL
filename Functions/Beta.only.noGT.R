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

#' Run Beta only with unknown rsnp GT, with inputs partially  made
#'
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. 
#' @param gene gene id for the gene to run
#' @param path.in path to input file
#' @param pattern pattern for input file
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param tag whether tags were used when inputs were prepared
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving files, if NULL gene_id.eqtl will be used
#' @param model path to stan model, defaults to '/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/beta.noGT.rsnp.priors.eff2.stan'
#' @keywords bayesian trecase unknown genotype regulatory snp
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from trecase model.
#' beta.only.nogt.with.in()

beta.only.nogt.with.in <- function(gene, path.in, pattern=".noGT.stan.in.rds", tag=c("no","yes"),covariates=1,out=".", prefix=NULL, model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/beta.noGT.rsnp.priors.eff2.stan') {

    
    if(covariates !=1){
        cov <- readRDS(covariates)
        if(!(is.matrix(cov))) stop("Covariates file is not a matrix")
       
    ## scale
        covariates <- apply(cov,2,scale,center=TRUE,scale=TRUE)
        rownames(covariates) <- rownames(cov)
        colnames(covariates) <- colnames(cov)
       
    }

    tag <- match.arg(tag)
    r.tag <- switch(tag=="yes","yes") ## to output results after stan, when tag.threshold is "no", returns NULL
    ## check if gene has input file

    in.file <- list.files(path.in, pattern=paste0(gene,pattern), full.names=T)

    if(length(in.file)==0){
        stop(cat("Gene ",gene, " has no input in ",path.in))
    } else {
      
        stan.noGT <- readRDS(in.file)

        ## test covariates and inputs have the same number of individuals
        samp <- sapply(stan.noGT, function(i) length(i$NB[[1]]))==nrow(covariates)
        
        if(!all(samp)){
            stop("Covariates file has different number of individuals compared to total counts in stan input")

            } else {
            
                stan.noGT2 <- lapply(stan.noGT, function(i) in.neg.beta.noGT.eff2(i, covar=covariates))     
    
                ## full model 
                mod <- stan_model(model)
                stan.full <-  mclapply(1:length(stan.noGT2),
                                       function (i) {
                                           s <- summary(sampling(mod,data=stan.noGT2[[i]], cores=1, refresh=0))$summary
                                           unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                                           return(s)
                                       })
                names(stan.full) <- names(stan.noGT)
                full.sum <- stan.bt(x=stan.full,y= "bj",rtag=r.tag,model="ase",gene=gene)
                if(!is.null(prefix)) {             
                    
                    write.table(full.sum, paste0(out,"/",prefix,".noGT.stan.summary.txt"), row.names=FALSE)
                } else {
                    
                    write.table(full.sum,paste0(out,"/",gene,".noGT.stan.summary.txt"), row.names=FALSE)
                }
                return(full.sum)  
            }
    }
    
}

