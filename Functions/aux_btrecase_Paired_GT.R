source('/home/ev250/Bayesian_inf/trecase/Functions/aux_btrecase_GT.R')

 
#' Aux for formatting 2 tretmenat inputs
#' @param l list with stan input for 2 treatments for 1 snps
#' @param dt data table to use for formatting
#' @param n name of elements on l to join
#' @export
#' @return data table with inputs for stan
#'
#' treat.form()

treat.form <- function(l, dt,n){
    tmp <- Reduce(function(a,b) merge(a,b, by=c("samples", "index"), all=T),
                  mapply( function(x,y) {
                      tmp2 <- data.table(samples=rep(dt$samples,dt[[x]]))
                      tmp2[, eval(n):=y[[n]]]
                      tmp2[ , index:=1:.N, by=samples]
                      return(tmp2)},
                      x=names(dt)[names(dt) !='samples'],
                      y=l,
                      SIMPLIFY=FALSE))
    ##change NA to 0 
               ### working here merging is incorrect, need an index                             
    tmp[is.na(tmp)] <- 0
    tmp[, index:=NULL]
    return(tmp)
}

#' Aux for returning genotype info
#' @param dt data table with inputs
#' @export
#' @return vector to return to stan
#' g.for()

g.for <- function(dt){
    n <- names(dt)[!names(dt) %in% 'samples']
    dt[, g:=get(n[1])][get(n[1])==0, g:=get(n[2])]
    return(dt$g)
}



    
#' Make inputs for running  Btrecase with known rsnp GT, optional refbias correction, preparing inputs more efficiently
#'
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to files with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R per treatment
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf files with ASE and GT for the chromosome where the gene is, for each treatment, same order as counts.f
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
#' btrecase.gt.paired.in()

btrecase.gt.paired.in <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,u.esnps=NULL, gene.coord,vcf,le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9, out=".", prefix=NULL, model=c("both","trecase","trec"), stan.model='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.refBias2.stan' ,
                                stan.trec='/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.eff2.stan',
                                prob=NULL, prior=NULL, ex.fsnp=NULL, AI_estimate=NULL, pretotalReads=100) {
    

   ## check inputs and extract inputs for gene
    model <- match.arg(model)
    if(!is.null(prior)){
        if(class(prior) != "list") stop("prior argument must be a list")
        if(any(!names(prior) %in% c("mean", "sd", "mix"))) stop("prior argument must be a named 'mean' and 'sd'")
        if(length(unique(sapply(prior, length))) !=1) stop("mean and sd for prior argument must have the same length")
    }
    
    ingene <- lapply(1:2, function(i) aux.in1(gene,
                                                   chr,
                                                   snps,
                                                   counts.f[i],
                                                   covariates,
                                                   e.snps,
                                                   u.esnps,
                                                   gene.coord,
                                                   vcf[i],
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
                                                   pretotalReads=pretotalReads))

    counts.g <- lapply(ingene, '[[', 'counts')
    covariates <- ingene[[1]]$covariates
    probs <- ingene[[1]]$probs

    rsnps.ex <- data.table(id=character(), reason=character())
    gcoord <- fread(gene.coord)

    ## get gt:ase counts
    gt.rs <- lapply(vcf, function(i) aux.snps(gene, chr, snps, gcoord, gene.coord, i,rsnps.ex))
    rs <- lapply(gt.rs, '[[', 'rs')
    gt.as <- copy(rs)
    rsnps.ex <- gt.rs[[1]]$rsnps.ex
    cis_window <- gt.rs[[1]]$cis

    ## QC vcf, make sure GT is the same in both files
    ## same names
    if(!Reduce(identical, lapply(gt.as, names))) stop("Names in vcf are not the same")
    ## check GT
    gt.cols <- names(gt.as[[1]])[-grep("_AS",names(gt.as[[1]]))]
    if(!Reduce(identical, lapply(gt.as, function(i) i[, gt.cols, with=F]))) stop("SNPs and/or genotypes are not the same in both vcf files")
    
    ## further process of rsnps
    rs <- rs[[1]][,grep("_AS",names(rs[[1]]),value=T):=NULL] # same for both treatments
    
    ## recode to 0,1,-1,2 scale
    rec.rs <- rec_mytrecase_rSNPs(x=rs$POS, y=rs)
    
    ## tagging
    if(tag.threshold!="no"){
        ## remove rsnps with var0, not possible to tag
        cols <- grep("_GT$", names(rec.rs),value=T)
        v <- apply(rec.rs[, cols, with=F],1,function(i)var(abs(i)))
        w <- which(v==0)
        if(length(w) == nrow(rec.rs)) stop("all SNPs have zero variance, cannot tag")
        id.tags <- help.tags(gene, tag.threshold, rs, prefix, rec.rs, out)
       
        if(length(w)){
            rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs[w, id], reason=rep("Snp with zero variance", length(w))))
        }
        rec.rs <- rec.rs[id %in% id.tags,]
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
            fsnps <- fsnps[id %in% gt.as[[1]]$id] # same in both treatments
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

                
            } else {

                ai <- NULL
            }
            
            
            if(nrow(fsnps)) {## stay
            
                if(is.character(snps)){ ## analysis on pre-specified set of rsnps
                    gt.as <- lapply(gt.as, function(i) i[id %in% c(rec.rs$id,fsnps$id),]) ## snps and fsnps only
                    ## make sure gt.as (here for fsnps) doesn't have missing values or unphased GT   
                    gt.qc <- lapply(gt.as, function(i) vcf.gt.qc(i))
                    s <- sapply(gt.qc, function(i) sum(i[c(2,4)])!=0)
                    if(any(s)) stop(cat("Invalid GT field for some exonic snps and samples \n",paste(names(gt.qc[[1]]), collapse=","), "\n",gt.qc[s], "\n","to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments \n", "if all exonic snps are removed the analysis will be done with total counts only"))           
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
                    w=which(rownames(rp)%in% gt.as[[1]]$id)
                    rp <- rp[w,, drop=FALSE]
                    ## make sure to select fsnps from reference panel
                    f.ase <- lapply(gt.as, function(i) i[id %in% fsnps$id,][id %in% rownames(rp),])                    
                    if(nrow(f.ase[[1]])==0) {
                        print("No fsnps ref panel")
                        rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs2$id, reason="No fsnps in reference panel"))
                    } else {
                        
                        ## get counts per hap for fsnp
                        counts.ase <- lapply(f.ase, tot.ase_counts)

                        if(!is.null(u.esnps)){
                            counts.ase <- lapply(counts.ase, function(i) aux.in2(gene, u.esnps, i, ex.fsnp, ai))
                            
                        }
                        
                        if(!all(sapply(counts.ase, function(i) is.character(i)))) { ## enough fsnps for ASE

                            ## check if sufficient indiviuals with ASE counts to proceed
                            s <- sapply(counts.ase, function(i) is.character(filt.fsnp(i,min.ase.snp=0)))
                            if(!all(s)) { ## Enough ind with ase counts
                                
                                print(paste("Effective number of fSNPs:", nrow(f.ase[[1]])))                        
                                
                                ## select rsnps in ref panel for full model
                                rs.full <- rs[id %in% rec.rs2$id, which(names(rs) %in% names(rec.rs2)), with=F]
                                
###################  run stan full model #######################
###### prepare stan inputs
                                ## reference panel fsnps and rsnps:
                                rp.f <- rp[f.ase[[1]]$id,,drop=FALSE]
                                rp.r <- rp[rs.full$id,,drop=FALSE]

                                ## get haplotype counts from fsnps

                                ## add AI_estimate per sample to stan inputs
                                if(!is.null(ai)) {
                                    ai <- ai[id %in% f.ase[[1]]$id,.(id, NREF_post, NALT_post,Total_post, AI_post)]
                                    ## get min AI_post
                                    min_AI <-  min(ai$AI_post)
                                }
                                

                                print("Preparing stan inputs")

                                ## order fsnps in couts.ase as in f.ase 

                                counts.ase=mapply(function(a,b) a[,unlist(lapply(b$id, function(i) grep(i, colnames(a), value = TRUE)))],
                                                  a=counts.ase,
                                                  b=f.ase,
                                                  SIMPLIFY=F)

                                stan.f <- mapply(function(a,b) stan.fsnp.noGT.eff(rp.f,a,b, NB="no", min.ase, min.ase.het, ai=ai),
                                                 a=f.ase,
                                                 b=counts.ase,
                                                 SIMPLIFY=F)

                                
                                if(!any(sapply(stan.f, is.character))) { ## proceed with full model
                                    
                                    counts <- lapply(counts.g, unlist) ## to avoid repeating in mclapply
                                    
                                    stan.in1<-mapply(function(a,b, c) mclapply(rs.full$id, function(i) stan.trecase.eff2(a, rp.1r=rp.r[i,,drop=FALSE], rp.f,b, rs.hap=rs.full[id==i,], rec.rsnp=rec.rs2[id==i,], c,min.ase, min.ase.het)),
                                                     a=counts,
                                                     b=f.ase,
                                                     c=stan.f,
                                                     SIMPLIFY=F)
                                    
                                    stan.in1 <- lapply(stan.in1, setNames, rs.full$id)
                                    

                                    ## remove rsnps with insufficient ase, min.ase.hets (haplotypes with fsnps not compatible with reference panel)
                                    w <- lapply(stan.in1, function(i) sapply(i, is.character))
                                    if(any(lapply(w,sum)>0)) {

                                        ## need to discard any rSNP that failed in any treatment
                                        failed <- Reduce(union, lapply(w, function(a) rs.full$id[a]))
                                                                       
                                        rsnps.ex <- rbind(rsnps.ex,data.table(id=failed, reason="Not enough het individuals with sufficient ase counts or Hap not in reference panel"))   
                                        
                                        stan.in1 <- lapply(stan.in1, function(i) i[!rs.full$id %in% failed])
                                    }
                                    if(length(stan.in1[[1]]) >= 1) {## same length both elements
                                        

                                        ## save stan.in1 to qc
                                        if(is.null(prefix)) {
                                            saveRDS(stan.in1, paste0(out,"/",gene,".GT.stan1.input.rds"))
                                            #saveRDS(colnames(counts.ase)[colSums(counts.ase)>5], paste0(out,"/",gene,".GT.fsnps.with.counts.rds"))
                                        } else{
                                            saveRDS(stan.in1, paste0(out,"/",prefix,".GT.stan1.input.rds"))
                                            #saveRDS(colnames(counts.ase)[colSums(counts.ase)>5], paste0(out,"/",prefix,".GT.fsnps.with.counts.rds"))
                                        }
                                        
                                        if(is.null(prior)){
                                            stan.in2 <- lapply(stan.in1, function(j) lapply(j, function(i) in.neg.beta.prob.eff2(i, covar=covariates)))
                                        } else {
                                            stan.in2 <- lapply(stan.in1, function(j) lapply(j, function(i) {
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
                                            ))
                                        }

                                        ## need to combine inputs from 2 treatments
                                        ## need to identify overlapping/distinct samples
                                        all <- Reduce(union, lapply(stan.f, function(i) names(i$m)))
                                        ## make indicator to link individuals with ASE counts
                                        ASEI <- do.call(data.table, lapply(stan.f, function(i) as.numeric(all %in%names( i$m))))
                                        ## add sample names to ASEI
                                        ASEI[, samples:=all]

                                        ## Combine inputs for each rsnp between the 2 treatments
                                        stan2t <- lapply(names(stan.in2[[1]]), function(i) {
                                            ##extract same snp per treatment
                                            l.sub <- lapply(stan.in2, '[[', i)
                                            ## ASE side                                            
                                            ## extract s and make data table to keep track of same individuals
                                            s.dt <- treat.form(l.sub, ASEI, "s")
                                            ## same with gase and pH
                                            gase.dt <- treat.form(l.sub, ASEI, 'gase')
                                            pH.dt <- treat.form(l.sub, s.dt, 'pH')
                                            ## prepare list with elements to return as inputs
                                            l <- mapply(function(l, dt,n) {
                                                tmp <- treat.form(l, dt, n)
                                                tmp[, samples:=NULL]
                                            },
                                                        dt=list(ASEI,s.dt),
                                                        n=c("m", "n"),
                                                        MoreArgs=list(l=l.sub),
                                            SIMPLIFY=FALSE)
                                            
                                            names(l) <- c("m", "n")
                                            if(!is.null(ai)){
                                                for(x in c('ai0', 'sdai0')){
                                                    l[[x]] <- treat.form(l.sub, s.dt, x)[, samples:=NULL]
                                                }
                                            }
                                            
                                            ## to return gase, pH and s
                                            for(x in c('gase', 'pH','s')){
                                                l[[x]] <- g.for(get(paste0(x, '.dt')))
                                            }
                                            l[['A']] <- length(all) ## union of ind with ASE counts across treatments
                                            l[['L']] <- length(l$pH)

                                            ## NB side and same elements in l.sub
                                            l[['Y']] <- Reduce(cbind, lapply(l.sub, function(i) i$Y))

                                            same <- c("N", "K", "g", "cov")
                                            if(!is.null(prior)) same <- c(same, 'k', 'aveP', 'sdP', 'mixP')
                                            for(x in same){
                                                l[[x]] <-l.sub[[1]][[x]]
                                            }
                                            return(l)
                                        })

                                        names(stan2t) <- names(stan.in1[[1]])
                                        
                                        ## count number of hets with sufficient ase counts to input in final report
                                        ASE.hets <- Reduce(function(a,b) paste(a,b, sep=","),
                                                           lapply(stan.in1, function(j) sapply(j, function(i) nrow(i$gm[abs(g.ase)==1,]))))
                                        ## get eaf for tags run in model
                                        eaf.t <- snp.eaf(le.file,names(stan2t),population)

                                        trecase.in <- list(trecase=stan2t, ASE.hets=ASE.hets, eaf.t=eaf.t, probs=probs, r.tag=r.tag, nhets=rec.rs2[id%in% names(stan2t),nhet],  nfsnps=nchar(unlist(strsplit(colnames(stan.f[[1]]$n[[1]]), ","))[1]))

                                        if(!is.null(AI_estimate)) trecase.in[['minAI']] <- min_AI
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
                
            }  else { ## closing from no AI estimates
                ## go to trec if model !=trecase, otherwise finish
                cat("No fSNPs with allelic imbalance estimates \n analysis will be done with total gene counts only")
                rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs$id, reason="No fSNPs with AI estimates"))
            }
            
        } ## closing from fsnps in gene
        if(model == "trecase") {
            if(nrow(rsnps.ex)){
                if(!is.null(prefix)) {             
                    write.table(rsnps.ex,paste0(out,"/",prefix,".trecase.excluded.rsnps.txt"), row.names=FALSE)
                    
                } else {
                    write.table(rsnps.ex,paste0(out,"/",gene,".trecase.excluded.rsnps.txt"), row.names=FALSE)
                }
            }
            if(exists("stan2t")){
                
                return(trecase.in)  
            } else {
                return("No rsnps were run with trecase")
            }
        }
        
    } else { ## when model==trec I add all rec.rs snps to prepare input

        rsnps.ex <- rbind(rsnps.ex,data.table(id=rec.rs$id, reason="requested model trec"))

    }
    
    ## neg.binom only 
    ## prepare input and run
    if(nrow(rsnps.ex)>0) {
        
        ex <- c("Missing or homo GT all samples", "zero variance", "rsnp with less than")
        w <- unlist(sapply(ex, function(i) grep(i,rsnps.ex$reason)))
        id.keep <- unique(rsnps.ex[!w,id])
        
        if(length(id.keep)) {
            in.neg <-lapply(counts.g, function(j)  mclapply(id.keep, function(i) input.neg.only.bj(DT1=j, DT2=rec.rs[id ==i,] ,covar=covariates)))

            in.neg <- lapply(in.neg, setNames, id.keep)

            neg2t <- lapply(id.keep, function(i) {
                ##extract  same snps
                l.sub <- lapply(in.neg, '[[', i)
                ## NB side and same elements in l.sub
                l <- list()
                l[['Y']] <- Reduce(cbind, lapply(l.sub, function(i) i$Y))

                same <- c("N", "K", "g", "cov")
                
                for(x in same){
                    l[[x]] <-l.sub[[1]][[x]]
                }
                
                if(!is.null(prior)) {
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
                }
                
                return(l)
            })
            names(neg2t) <- id.keep
            
            ## get eaf for tags run in model
            eaf.t.neg <- snp.eaf(le.file,names(neg2t),population)
            ## if tag not in ref panel (possible with neg binom side)
            if(length(id.keep) > nrow(eaf.t.neg)){
                eaf.t.neg <- rbind(eaf.t.neg, data.table(snp=id.keep[which(!id.keep %in%  eaf.t.neg$snp)], eaf=NA))
                ## sort as in id.keep
                eaf.t.neg <- eaf.t.neg[order(match(snp, id.keep))]
            }
            if(!is.null(prefix)){             
                write.table(rsnps.ex,paste0(out,"/",prefix,".trecase.excluded.rsnps.txt"), row.names=FALSE)
                
                
            } else {
                write.table(rsnps.ex,paste0(out,"/",gene,".eqtl.trecase.excluded.rsnps.txt"), row.names=FALSE)
                
                
            }
            neg.in <- list(neg=neg2t, eaf.t=eaf.t.neg, r.tag=r.tag, probs=probs,  nhets=rec.rs[id %in% id.keep,nhet])
            if(exists("stan2t")){
                return(list(trecase=trecase.in, neg=neg.in))
            } else{
                return(list(neg=neg.in))
            }
            

        }
    } else {
        if(model == "trec") {
            return("No rsnps can be run with trec")
        }

        if(exists("stan2t")){
            if(!is.null(prefix)) {             
                
                write.table(rsnps.ex,paste0(out,"/",prefix,".trecase.excluded.rsnps.txt"), row.names=FALSE)
            } else {
                
                write.table(rsnps.ex,paste0(out,"/",gene,".eqtl.trecase.excluded.rsnps.txt"), row.names=FALSE)
                
            }
            return(trecase=trecase.in)
            
        } else {
            return("No snps  can be run with trecase or trec")
        }
    }
    if(!exists("stan2t")){
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

