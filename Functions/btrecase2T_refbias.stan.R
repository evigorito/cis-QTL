##23/02/2020, Improves from Btrecse2T.R by allowing choice of checking het.f by tissue or joint and by jointly testing info score across all samples.
suppressMessages(library(ggplot2))
suppressMessages(library(MASS))
library(parallel)
suppressMessages(library(rstan))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


source('/home/ev250/Bayesian_inf/trecase/Functions/aux_btrecase2T.R')

library(GUESSFM)

#' Run Btrecase for 2 tissues in one model but checking heterogeneity independently.
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. 
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the exonic snps for the chromosome where the gene is
#'  @param sample.file sample file for the reference panel (sample description), to be used if ex.fsnp test is required and population is not the whole reference panel
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL, defaults to EUR
#' @param maf cut-off for maf, defaults to 0.05
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.snp minum number of ASE counts for a single snp to be considered, for a particular individual
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving files, if NULL gene_id.eqtl will be used
#' @param model path to stan model, defaults to /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Scripts/ba01bd01refbias.stan
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.95 and 0.99 quantiles
#' @param prior named list: mean= vector with the mean of Gaussians, sd= vector with Gaussians sd for eQTL effect prior, mix=vector with mixing proportions. Deffaults to NULL
#' @param ex.fsnp, p-value cut-off for fSNPs to exclude, defaults to NULL
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, defaults to NULL
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to NULL
#' @param skin character vector with the 2 Tissues to study
#' @param fishjoint whether to run Fisher test jointly in all samples or by skin, for QC purposes for comparing with 1T model defaults to by skin
#' @keywords bayesian trecase unknown genotype regulatory snp
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from trecase model.
#' btrecase2T.stan()

btrecase2T.stan <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,u.esnps=NULL,gene.coord,vcf,sample.file=NULL, le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), maf=0.05, min.ase=5,min.ase.snp=5,min.ase.n=5,tag.threshold=.9, info=0.3, out=".", prefix=NULL, model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/psoriasis/run_Btrecase/Scripts/2TrefbiasV2SpikeMix.stan', prob=NULL, prior=NULL, ex.fsnp=NULL, AI_estimate=NULL,pretotalReads=NULL, skin, fishjoin=NULL) {

    if(!file.exists(model)) stop("Invalid stan model name")
    inputs <- btrecase.nogt.rna.refbias2T.In(gene,
                                             chr,
                                             snps,
                                             counts.f,
                                             covariates,
                                             e.snps,
                                             u.esnps,
                                             gene.coord,
                                             vcf,
                                             sample.file,
                                             le.file,
                                             h.file,
                                             population,
                                             maf,
                                             min.ase,
                                             min.ase.snp,
                                             min.ase.n,
                                             tag.threshold,
                                             info,
                                             out,
                                             prefix,
                                             ex.fsnp,
                                             prob,
                                             AI_estimate,
                                             pretotalReads,
                                             skin,
                                             fishjoin)

    if(is.character(inputs)) stop(inputs)
    ## rename inputs
    probs <- Reduce(intersect, inputs$probs)
    r.tag <- inputs$r.tag
    het.f <- inputs$het.f
    nfsnps <- inputs$nfsnps
    info.ok <- inputs$inp$info.ok
    stan2T <- inputs$inp$stan.in
    if(!is.null(AI_estimate)) min_AI <- inputs$inp$min_AI
    
    
    cat("Running stan for ", length(info.ok), " rSNPS with model\n", model)
    mod <- stan_model(model)
    
     stan.full <-  mclapply(1:length(stan2T),
                           function (i) {
                               post <- sampling(mod,data=stan2T[[i]], cores=1, refresh=0)
                               s <- summary(post, pars=c("ba", "bd", "bp", "bn"), use_cache=F, probs=probs)$summary[1:4,]
                               
                               unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                              
                               ## process summary
                               dt <- data.table(s, keep.rownames=T)
                               quant <- grep("%", names(dt), value=T)
                               setnames(dt, "rn", "Param")
                               ## add Signif column
                               dt[, Signif:="no"][sign(get(quant[1])) == sign(get(quant[length(quant)])) , Signif:="yes"]
                               dt[, rSNP:=names(stan2T)[i]]

                                ## extract params from posterior to get PEP
                               e <- rstan::extract(post, pars=c("ba", "bd","bp", "bn"))
                               post.neg <- sapply(e, function(i) sum(i<0)/length(i))
                               ## PEP: posterior area opposite sign to mean posterior
                               dt[, PEP:=post.neg][mean<0, PEP:= 1-PEP]
                               
                               dt.wide <- reshape(dt, idvar="rSNP", timevar="Param", direction="wide")
                               ##print(dt.wide)
                               ##return(list(dt=dt.wide, edt=edt))
                               return(dt.wide)
                           })

    ## get maf for tags run in model
    maf.t <- snp.eaf(le.file,names(stan2T),population)       
    
    full.sum <- stan.2T(x=stan.full, rtag=r.tag ,gene=gene,EAF=maf.t, info=info.ok,
                        nfsnps=nfsnps,
                        min.pval= het.f,
                        probs=probs)

    if(!is.null(AI_estimate)) full.sum[, min_AI:= min_AI]

    if(!is.null(prefix)){
        write.table(dt,paste0(out,"/",prefix,".2Tissues.summary.txt"), row.names=FALSE)
    } else {
       write.table(full.sum,paste0(out,"/",gene,".2Tissues.summary.txt"), row.names=FALSE)

    }
}
