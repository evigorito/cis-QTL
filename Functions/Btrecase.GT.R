
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

#' Run Btrecase with known rsnp GT 
#'
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the chromosome where the gene is
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param nhets minimun number of het individuals in order to run the minumn model (NB only), defaults to 5
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.het minimun number of het individuals with the minimun of ASE counts in order to run the ASE side of the model, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param q.test whether to select snps to run based on a pre-filtering test, character string indicating "yes", "no" (default)
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving tables, if NULL gene_id.eqtl will be used
#' @param model whether to run trecase, trec or both, defaults to both
#' @keywords bayesian trecase known genotype regulatory snp
#' @export
#' @return data.table with summary of gene-snp associations. Saves the summary table in "out" dir as /out/prefix.main.txt. When using tags, saves /out/prefix.tags.lookup.txt. Saves a table of excluded rsnps from trecase model.
#' btrecase.gt()

btrecase.gt <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,gene.coord,vcf,le.file,h.file,nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no", out=".", prefix=NULL, model=c("both","trecase","trec")) {
    
    ## check inputs
    model <- match.arg(model)
    files <- c(counts.f,e.snps,gene.coord,vcf,le.file,h.file)
    w <- !file.exists(files)
    if(any(w)) stop(paste("invalid file names ", paste(files[w], collapse= ", ")))
    if(!(tag.threshold >=0 & tag.threshold <=1 | tag.threshold =="no")) stop("invalid  'tag.threshold' argument")
    na.q.test <- pmatch(q.test, c("yes","no"))
    if(is.na(na.q.test)) stop("invalid 'q.method' argument")
    if(!(dir.exists(out))) stop("invalid 'out' argument")
    num <- list(nhets,min.ase,min.ase.het)
    names(num) <-  c("nhets","min.ase","min.ase.het")
    w <- !sapply(num, is.numeric)
    if(any(w)) stop(paste("invalid arguments:",paste(names(num)[w], collapse=", ")))
    if(!(covariates==1 | file.exists(covariates))) stop("invalid path for covariates file")
    
    ## Extract inputs for gene

    ## get counts 
    counts.g <- fread(paste("grep -e gene_id -e ",gene,counts.f), header=TRUE)
    if(nrow(counts.g)==0) stop("Gene id is not found in count matrix")
    ## get covariates and scale
    if(covariates !=1){
        covariates <- readRDS(covariates)
        if(!(is.matrix(covariates))) stop("Covariates file is not a matrix")
        if(nrow(covariates)!=ncol(counts.g)-1)  stop(cat("Number of individuals in covariates matrix:", nrow(covariates), "\n Number of individuals in gene counts", ncol(counts.g)-1, "\n please adjust"))
        ## scale
        covariates <- apply(covariates,2,scale,center=TRUE,scale=TRUE)
    }
    
    ## get rsnps and extract GT and ASE, remove non-informative snps (missing or hom in all samples)
    ## create dt to collect rsnps excluded from analysis
    rsnps.ex <- data.table(id=character(), reason=character())
    if(is.numeric(snps)){
        cis_window <- cl_coord(gene.coord,22,gene,cw=snps)
        gt.as <- vcf_w(vcf,22, st=cis_window["start"], end=cis_window["end"], exclude="yes")
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
        st_end <- cl_coord(gene.coord,chr,gene,cw=0)
        ## construct cis-window with snps, making sure to include the whole gene
        m <- min(pos) < st_end[1]
        M <- max(pos) > st_end[2]
        cis_window <- ifelse(m, min(pos), st_end[1])
        cis_window <- c(cis_window, ifelse(M, max(pos), st_end[2]))
        gt.as <- vcf_w(vcf,22, cis_window["start"], cis_window["end"], exclude = "yes")
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
    ## recode to 0,1,-1,2 scale
    rec.rs <- rec_mytrecase_rSNPs(x=rs$POS, y=rs)
    if(tag.threshold!="no") {
        ## Group rsnps by r2, recode rec.rs for input in tags function from GUESSFM
        if(nrow(rs)==1) stop("Only one regulatory snp to test, please set tag.threshold='no' \n Cannot cluster one snp only")
        re.guess <- rec.guess(DT=rec.rs)
        x <- as(re.guess-1, "SnpMatrix")
        rtag <- tag(X=x,tag.threshold=tag.threshold)
        ## save rtag as data.table 
        dt <- data.table(tag=tags(rtag), SNP=rtag@.Data)
        if(!is.null(prefix)){
            write.table(dt,paste0(out,"/",prefix,".tags.lookup.txt"), row.names=FALSE)
        } else {
            write.table(dt,paste0(out,"/",gene,"eqtl.tags.lookup.txt"), row.names=FALSE)
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
    fsnps <- tryCatch({fread(paste("grep", gene, e.snps))}, error=function(e) {paste("No entry for gene",gene, "in",e.snps)})
    
    if(isTRUE(grep("No entry",fsnps)==1)) { ## no fsnps
        cat(fsnps, "\n Analysis will be done with total gene counts only")
        rsnsp.ex <- rbind(rsnps.ex, data.table(id=rec.rs$id, reason=fsnps))
    } else {     
        fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
        fsnps <- fsnps[V3 %in% gt.as$ID]
        if(is.character(snps)){ ## analysis on pre-specified set of rsnps
            gt.as <- gt.as[id %in% c(rec.rs$id,fsnps$id),] ## snps and fsnps only
            ## make sure gt.as (here for fsnps) doesn't have missing values or unphased GT   
            gt.qc <- vcf.gt.qc(gt.as)
            if(sum(gt.qc[c(2,4)])!=0) stop(cat("Invalid GT field for some exonic snps and samples \n",paste(names(gt.qc), collapse=","), "\n",gt.qc, "\n","to remove snps or samples with wrong GT use vcf.gt.qc with appropiate arguments \n", "if all exonic snps are removed the analysis will be done with total counts only"))           
        }
        
        ## get info from reference panel 
        ## matrix for reference panel haps
        rp <- haps.range(file1=le.file,file2=h.file,cis_window)
        ## keep record of rsnps not in reference panel
        w <- which(!rec.rs$id %in% rownames(rp))
        rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs$id[w],reason=rep("Not in reference panel", length(w))))
        rec.rs2 <- rec.rs[!w,] ## for full model, only rsnps in ref panel
        w <- which(rownames(rp)%in%  rec.rs2$id)
        if(length(w)!=0) {  ## proceed with full model, neg only for rsnps not in ref panel           
            w=which(rownames(rp)%in% gt.as$id)
            rp <- rp[w,, drop=FALSE]
            ## make sure to select fsnps from reference panel
            f.ase <- gt.as[ID %in% fsnps$V3,]
            f.ase[, id:= paste(POS,REF,ALT, sep=":")]
            f.ase <- f.ase[id %in% rownames(rp),]
            if(nrow(f.ase)==0){
                print("No fsnps ref panel")
                rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs2$id, reason="No fsnps in reference panel"))
            } else {
                print(paste("Effective number of fSNPs:", nrow(f.ase)))
                ## get counts per hap for fsnp
                counts.ase <- tot.ase_counts(x=f.ase, y=counts.g , z=rec.rs2)
                names(counts.ase) <- rec.rs2$id
                ## select rsnps with at least min.ase.het  with min.ase counts
                filt.counts <- lapply(counts.ase, filt.rsnp,ase=min.ase,n=min.ase.het)
                names(filt.counts) <- rec.rs2$id
                ## select data tables from filt.counts
                w <- sapply(filt.counts, is.data.table)
                if(!all(w)) { 
                    rsnps.ex <- rbind(rsnps.ex,data.table(id=rec.rs2$id[!w], reason=paste("less than",min.ase.het,"het ind. with sufficient ase counts")))
                    filt.counts <- filt.counts[w]    
                }
                if(length(filt.counts)>0){ ## full model
                    
                    ## Apply rule to eliminate null associations
                    if(q.test=="yes"){
                        q.test <- lapply(filt.counts, q.rsnp,ase=min.ase,n=min.ase.het)
                        q.test <- q.test[which(sapply(q.test, is.numeric) == TRUE)]
                        keep <- lapply(q.test, function(i) abs(i[[1]])>= 0.2 | abs(i[[2]]) >= 2)      
                        keep <- keep[keep==TRUE]            
                        ## remove snps from filt.counts
                        w <- which(names(filt.counts) %in% names(keep))
                        if(length(w)<length(filt.counts)){
                            rsnps.ex <- rbind(rsnps.ex, data.table(id=names(filt.counts)[-w], reason="below cut-off for qtest"))
                            filt.counts <- filt.counts[w]
                        }
                    }
                    
                    ## rsnps for full model
                    rec.full <- rec.rs2[id %in% names(filt.counts)]
                    rs.full <- rs[id %in% rec.full$id, which(names(rs) %in% names(rec.full)), with=F]
                    
###################  run stan full model #######################
###### prepare stan inputs
                    ## reference panel fsnps and rsnps:
                    rp.f <- rp[f.ase$id,,drop=FALSE]
                    rp.r <- rp[rs.full$id,,drop=FALSE]
                    
                    rp2 <- lapply(1:nrow(rp.r), function(i) t(rbind(rp.f,rp.r[i,]))) ## get ref panel fsnps and rsnp in the same format I have functions from simulations
                    names(rp2) <- rs.full$id
                    
                    ## calculate P(H|G) ref panel for each rsnp                   
                    rp.hap.pairs <- mclapply(rp2, p.hap.pair) 

                    ## get haps for fsnps and each rsnp
                    h.samp <- mclapply(rec.full$id, function(i) hap_sam(x=f.ase,y=i,z=rs.full))

                    ## select filt.counts with snps in ref panel
                    filt.counts <- filt.counts[which(names(filt.counts) %in% rec.full$id)]
                    
                    stan.in1 <- mclapply(1:length(h.samp), function(i) stan.neg.beta.prob.eff(g=h.samp[[i]][[1]] + h.samp[[i]][[2]] , p.hap.pairs=rp.hap.pairs[[i]], h1=h.samp[[i]][[1]], h2=h.samp[[i]][[2]], geno.exp=filt.counts[[i]], ase=min.ase, n=min.ase.het ))
                    names(stan.in1) <- rec.full$id

                    ## remove rsnps with not sufficient ase, min.ase.hets (haplotypes with fsnps not compatible with reference panel)
                    w <- which(stan.in1 == "Not enough individuals with ASE counts")
                    if(length(w)!=0){
                        rsnps.ex <- rbind(rsnps.ex,data.table(id=rec.full$id[w], reason="Not enough individuals ase counts, Hap not in reference panel"))   
                        
                        stan.in1 <- stan.in1[-w]
                    }
                        if(length(stan.in1) >= 1){
                            
                            stan.in2 <- lapply(stan.in1, function(i) in.neg.beta.prob.eff2(i, covar=covariates))

                            mod <- stan_model('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.beta.prob.phasing.priors.eff2.stan')
                            stan.full <-  mclapply(1:length(stan.in2),
                                                   function (i) {
                                                     
                                                        s <- summary(sampling(mod,data=stan.in2[[i]], cores=1, refresh=0))$summary
                                                       unload.ddl(mod) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                                                       return(s)
                                                   })
                            names(stan.full) <- names(stan.in1)
                            ## count number of hets with sufficient ase counts to input in final report
                            ASE.hets <- sapply(stan.in1, function(i) nrow(i$gm[abs(g.ase)==1,]))
                            full.sum <- stan.bt(x=stan.full,y= "bj",rtag=r.tag,model="trec-ase", nhets=rec.full[id%in% names(stan.full),nhet],ASE.het=ASE.hets)
                            ##cat("number of rsnp full model: ",length(stan.full))
                        
                        } ## closing from no fnsps ref panel
                    } ## closing from no rnsps ref panel
                }
            }
        }
    
        if(model == "trecase"){
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
                } else {
                    
                    write.table(full.sum,paste0(out,"/",gene,".stan.summary.txt"), row.names=FALSE)
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
        id.keep <- rsnps.ex[!w,id]
        if(!!length(id.keep)){            
            in.neg <- mclapply(id.keep, function(i) input.neg.only.bj(DT1=counts.g[,2:ncol(counts.g),with=FALSE], DT2=rec.rs[id ==i,] ,covar=covariates))
            
            mod2 <- stan_model('/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/neg.only.eff2.stan')
            stan.neg <-  mclapply(in.neg, function (i) {
                s <- summary(sampling(mod2,data=i, cores=1, refresh=0))$summary
                unload.ddl(mod2) ##removing unnecessary dlls, when they reach 100 R gives error https://github.com/stan-dev/rstan/issues/448
                return(s)
            })
            names(stan.neg) <- id.keep         
            neg.sum <- stan.bt(x=stan.neg, y="bj",rtag=r.tag,model="trec", nhets=rec.rs[id %in% id.keep,nhet])
            if( exists("full.sum")){
                neg.sum <- rbind(full.sum, neg.sum)
                setorder(neg.sum,`log(aFC)_null`,-`log(aFC)_d.aux`, -model)
                
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
        if(!exists("full.sum")){
                    if(!is.null(prefix)) {             
                        write.table(full.sum, paste0(out,"/",prefix,".stan.summary.txt"), row.names=FALSE)
                    } else {
                        write.table(full.sum,paste0(out,"/",gene,".stan.summary.txt"), row.names=FALSE)
                    }
                    return(full.sum)  
                } else {         
                    return("No snps run using trecase or trec")
                }    
    }
}










