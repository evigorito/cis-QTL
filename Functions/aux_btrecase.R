source('/home/ev250/Cincinatti/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')

source('/home/ev250/Bayesian_inf/trecase/Functions/stan.eff.R')

##library(devtools)

##install_github("chr1swallace/GUESSFM", ref="groups")
library(GUESSFM)


#' First aux function to prepare inputs to run Btrecase with or without known rsnp GT: check inputs
#'
#' This function allows you to check arguments with and without GT
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param covariates path to matrix of covariates prepared in inputs.R, if no covariates, covariates =1, default
#' @param e.snps path to file listing exonic snps for the chromosome where the gene is, prepared in input.R
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with ASE and GT for the chromosome where the gene is
#' @param sample.file sample file for the chromosome of interest for the reference panel (sample description), to be used if ex.fsnp test is required and populatios is not the whole reference panel
#' @param le.file path to gz legend file (legend.gz) for the chromosome of interest for the reference panel (snp description)
#' @param h.file path to gz haplotpe file for the chromosome of interest for the reference panel (haplotypes for all samples in reference panel)
#' @param population ethnicity to get EAF for rsnp: AFR AMR EAS EUR SAS ALL, defaults to EUR. It is also used to exclude fSNPs if frequency is different from reference panel when using no GT and selecting ex.fsnps argument, except if population=ALL
#' @param maf cut-off for maf, defaults to NULL as only used with noGT
#' @param nhets minimun number of het individuals in order to run the minumn model (NB only), defaults to 5
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.het minimun number of het individuals with the minimun of ASE counts in order to run the ASE side of the model, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving tables, if NULL gene_id.eqtl will be used
#' @param model GT: whether to run trecase, trec or both, noGT NULL, defaults to NULL
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.99 and 0.95 quantiles
#' @param ex.fsnp, if character: vector with pos:ref:alt for fsnps to exclude,  defaults to NULL
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, defaults to NULL
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @keywords bayesian trecase known genotype regulatory snp
#' @export
#' @return list with counts, covariates and probs argument to prepare inputs
#' aux.in1()

aux.in1 <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,u.esnps=NULL, gene.coord,vcf,sample.file=NULL, le.file, h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), maf=NULL, nhets=NULL,min.ase=5,min.ase.het=NULL, min.ase.snp=NULL, min.ase.n=NULL,tag.threshold=.9, info=NULL, out=".", model=NULL, stan.model=NULL, prob=NULL, ex.fsnp=NULL, AI_estimate=NULL, pretotalReads=100) {

    ## check inputs:
    if(!is.null(model)) {
        w <- model %in% c("both","trecase","trec")
        if (!all(w)) stop("Not valid model selection")
    }
    
    if(!(dir.exists(out))) stop("invalid 'out' argument")
    population <- match.arg(population)
    files <- c(counts.f,e.snps,gene.coord,vcf,le.file,h.file)
    if(!is.null(u.esnps)) files <- c(files, u.esnps)
    if(!is.null(sample.file)) files <- c(files,sample.file)
    if(!is.null(stan.model)) files <- c(files, stan.model)
    if(!is.null(AI_estimate)){
        files <- c(files,AI_estimate)
        if(!is.numeric(pretotalReads)) stop("pretotalReads requires a numeric value")
    }
     w <- !file.exists(files)
    if(any(w)) stop(paste("invalid file names ", paste(files[w], collapse= ", ")))

    if(!(tag.threshold >=0 & tag.threshold <=1 | tag.threshold =="no")) stop("invalid  'tag.threshold' argument")
    
    num <- list(min.ase=min.ase)

    if(!is.null(maf)) num <- c(num, maf=maf)
    if(!is.null(nhets)) num <- c(num, nhets=nhets)
    if(!is.null(min.ase.het)) num <- c(num, min.ase.het)
    if(!is.null(min.ase.snp)) num <- c(num, min.ase.snp=min.ase.snp)
    if(!is.null(min.ase.n)) num <-  c(num, min.ase.n=min.ase.n)
    if(!is.null(info)) num  <- c(num, info=info)
    w <- !sapply(num, is.numeric)
    if(any(w)) stop(paste("invalid arguments:",paste(names(num)[w], collapse=", ")))
    
    if(covariates!=1) {
       if(!file.exists(covariates)) stop("invalid names for covariates file")
    }
    if(is.character(snps)){
        pos <- as.numeric(sapply(strsplit(snps, split=":"), `[[`,1))
        w <- which(!is.na(pos))
        if(!length(w)) stop(cat("Invalid format for snps ", snps[w]))
    }

    ## set default quantile for stan posterior to get 99 and 95 CI
    probs <- if(is.null(prob)) c(0.005, 0.025, 0.25, 0.50, 0.75, 0.975, 0.995) else c((1-prob)/2, 0.25, 0.5, 0.75, (1+prob)/2)

    ## Extract inputs for gene

    ## get counts 
    counts.g <- fread(cmd=paste("grep -e gene_id -e ",gene,counts.f), header=TRUE)
    if(nrow(counts.g)==0) stop("Gene id is not found in count matrix")
    counts.g <- counts.g[,2:ncol(counts.g),with=F] ## removes gene_id
    ## get covariates and scale
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
    
    
    return(list(counts=counts.g, covariates=covariates, probs=probs))
    

}


#' Second aux function to prepare inputs to run Btrecase with or without known rsnp GT: ASE counts mapping unique fSNPs
#'
#' This function allows you to check arguments with and without GT
#' @param gene gene id for the gene to run
#' @param u.esnps whether to use unique exonic snps per gene, defaults to NULL when it is not necessary if strand info is known
#' @param case, data table with ASE counts
#' @param ex.fsnp character vector with id of fsnps to exclude
#' @param ai data table with allelic imbalance estimates
#' @keywords bayesian trecase aux input function
#' @export
#' @return data table with filtered ASE counts
#' aux.in2()

aux.in2 <- function(gene, u.esnps,case,ex.fsnp=NULL, ai=NULL) {
    
     ufsnps <- tryCatch({fread(cmd=paste("grep", gene, u.esnps))}, error=function(e) {paste("No entry for gene",gene, "in",u.esnps)})
                if(!(isTRUE(grep("No entry",ufsnps)==1) | nrow(ufsnps) == 0) ) { ## ASE check
                    
                    ufsnps[ ,id:= paste(V2, V4,V5, sep=":")]
                    ## exclude fsnps if required

                    if(!is.null(ex.fsnp)){
                        ufsnps <- ufsnps[!id %in% ex.fsnp,]
                    }

                    if(!is.null(ai)) ufsnps <- ufsnps[ id %in% ai$id, ]

                    ## check if any unique fSNP to use
                    if(!nrow(ufsnps)) {
                        case=paste("No unique fSNPs for  gene",gene)
                    }
                    
                    keep <- unlist(lapply(ufsnps$id, function(i) grep(i, colnames(case), value=T)))

                    ## check if anything to keep
                    if(!length(keep)){
                        case=paste("No unique fSNPs for  gene",gene)
                    }
                    
                    rem <- colnames(case)[which(!colnames(case) %in% keep)]
                    ## make counts to 0, I dont remove the fSNPs because I want to use them for phasing
                    if(length(rem)) case[,rem]  <-  0
                   
                } else {
                    ## make case character so flow goes to NB only
                    case= paste("No unique fSNPs for  gene",gene)
                    
                }
    return(case)

}


#' Third aux function to prepare inputs to run Btrecase without known rsnp GT: prepare stan inputs wrap
#'
#' This function allows you to prepare stan inputs final stages
#' @param gene gene id for the gene to run
#' @param ai estimates for allelic imbalance, defaults to NULL
#' @param case, data table with ASE counts
#' @param rp.f reference panel haps for fSNPs
#' @param rp.r ref panel haps for rSNPs
#' @param f.ase data table with genotypes for fSNPs
#' @param counts.g data table with total counts for gene
#' @param covariates 
#' @param min.ase
#' @param min.ase.n
#' @param info score
#' @param snps.ex data table with rsnps to exclude from running
#' @param prefix character with prefix to add for saving files, defaults to NULL
#' @param out path to save outputs, default to current directory
#' @keywords bayesian trecase aux input function
#' @export
#' @return list with stan input
#' aux.in3()

aux.in3 <- function(gene, ai=NULL, case, rp.f, rp.r, f.ase, counts.g, covariates, min.ase=5, min.ase.n=5, info=0.3, snps.ex, prefix=NULL, out='.') {
    
    stan.f <- fsnp.prep2(rp.f, f.ase, case , min.ase, min.ase.n, ai)

    
    if(is.character(stan.f)) return(stan.f)
                       
    stan.noGT <- mclapply(1:nrow(rp.r), function(i) stan.trecase.rna.noGT.eff2(counts.g, rp.1r=rp.r[i,,drop=FALSE], rp.f, stan.f))

    names(stan.noGT) <- rownames(rp.r)
    
    ## remove rsnps when genotypes of fsnps are not compatible with reference panel
    w <- sapply(stan.noGT, is.character)

    if(any(w)){
        snps.ex <- rbind(snps.ex,data.table(id=names(stan.noGT)[w], reason=unlist(stan.noGT[w])))   
        stan.noGT <- stan.noGT[!w]
        rp.r <- rp.r[names(stan.noGT),]
    }
    
    if(all(w)) return("Genotypes of fSNPs are not compatible with reference panel")
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

    if(length(stan.noGT) >= 1) {

        ## save inputs for QC
        if(is.null(prefix)) {
            saveRDS(stan.noGT, paste0(out,"/",gene,".noGT.stan.input.rds"))
        } else{
            saveRDS(stan.noGT, paste0(out,"/",prefix,".noGT.stan.input.rds"))
        }
        
        stan.noGT2 <- lapply(stan.noGT, function(i) {
            l=in.neg.beta.noGT.eff2(i, covar=covariates[names(i$NB$counts),, drop=F])
            return(l)
        })

        ## remove het.f from here, use the calculation in btrecase.nogt.rna.refbias.In
        ##het.f <- suppressWarnings(prop_het(f.ase,rp.f,gene)) ## missing values converted to NA
        nfsnps=length(unique(Reduce(c, stan.f$f.comb)))

        ## return(list(stan.noGT2=stan.noGT2, het.f=het.f, nfsnps=nfsnps, info.ok=info.ok, snps.ex=snps.ex))
        return(list(stan.noGT2=stan.noGT2, nfsnps=nfsnps, info.ok=info.ok, snps.ex=snps.ex))

    } else {

        return("None of the snps  met the conditions to be run by model")
    }

}

    
#' Prepare inputs for Btrecase with unknown rsnp GT and missing values for GT fsnps with reference panel bias correction
#'
#' This function allows you to get inputs for Btrecase for one gene and multiple pre-selected snps. 
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
#' @param population ethnicity to set a cut-off for maf: AFR AMR EAS EUR SAS ALL, defaults to EUR
#' @param maf cut-off for maf, defaults to 0.05
#' @param min.ase minimun number of ASE counts for an individual in order to be included, defaults to 5
#' @param min.ase.snp minum number of ASE counts for a single snp to be considered, for a particular individual
#' @param min.ase.n minimun number individuals with the minimun of ASE counts, defaults to 5.
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param info numeric cut-off for var(E(G))/var(G), var(E(G)) is the expected variance for input and var(G) for reference panel, similar to info score in impute2, defaults to 0.3. rsnps with lower info wont be run by stan.
#' @param out path to save outputs, default to current directory
#' @param model GT: whether to run trecase, trec or both, noGT NULL, defaults to NULL
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.95
#' @param ex.fsnp, if character: vector with pos:ref:alt for fsnps to exclude, if numeric  p-value cut-off for fSNPs to exclude,  defaults to NULL
#' @param prob  number p∈(0,1) indicating the desired probability mass to include in the intervals, defaults to 0.95
#' @param AI_estimate full name to data table with AI estimates for reference panel bias for fSNPs, defaults to NULL
#' @param pretotalReads numeric indicating a cut-off for total initial reads to consider AI estimates, defaults to 100
#' @keywords bayesian trecase unknown genotype regulatory snp reference panel bias
#' @export
#' @return list with 1)c.ase and 2)stan.noGT object
#' btrecase.nogt.rna.refbias.In()

btrecase.nogt.rna.refbias.In <- function(gene, chr, snps=5*10^5,counts.f,covariates=1,e.snps,u.esnps=NULL,gene.coord,vcf,sample.file=NULL, le.file,h.file,population=c("EUR","AFR", "AMR", "EAS",  "SAS", "ALL"), maf=0.05, min.ase=5,min.ase.snp=5,min.ase.n=5,tag.threshold=.9, info=0.3, out=".", model =NULL, prefix=NULL, ex.fsnp=NULL, prob=NULL, AI_estimate=NULL, pretotalReads=100) {
  
    ## check inputs and extract inputs for gene
    
    ingene <- aux.in1(gene,
                      chr,
                      snps=snps,
                      counts.f,
                      covariates=covariates,
                      e.snps,
                      u.esnps,
                      gene.coord,
                      vcf,
                      sample.file=sample.file,
                      le.file,
                      h.file,
                      population=population,
                      maf=maf,
                      nhets=NULL,
                      min.ase=min.ase,
                      min.ase.het=NULL,
                      min.ase.snp=min.ase.snp,
                      min.ase.n=min.ase.n,
                      tag.threshold=tag.threshold,
                      info=info,
                      out=out,
                      model=NULL,
                      prob=prob,
                      ex.fsnp=ex.fsnp,
                      AI_estimate=AI_estimate,
                      pretotalReads=pretotalReads)
 
   
    counts.g <- ingene$counts
    covariates <- ingene$covariates
    probs <- ingene$probs

    
    ## create dt to collect snps excluded from analysis
    snps.ex <- data.table(id=character(), reason=character())
    
    ## get fSNPs (feature snps or exonic snps)
    fsnps <- tryCatch({fread(cmd=paste("grep", gene, e.snps))}, error=function(e) {paste("No entry for gene",gene, "in",e.snps)})


    
    if(isTRUE(grep("No entry",fsnps)==1) | nrow(fsnps) == 0) { ## no fsnps
        stop(paste("No entry for gene",gene, "in",e.snps))
    }

    fsnps[ ,id:= paste(V2, V4,V5, sep=":")]
    ## exclude fsnps if required

    if(!is.null(ex.fsnp)){
        fsnps <- fsnps[!id %in% ex.fsnp,]
    }

    if(nrow(fsnps)==0)  stop(paste("No valid fSNPs after excluding ", ex.fsnp))
   
    ## get fsnps and extract GT and ASE, remove non-informative snps (missing or hom in all samples)
    fsnps.window <- c(min(fsnps$V2),max(fsnps$V2))
    gt.as <- vcf_w(vcf,chr, st=fsnps.window[1], end=fsnps.window[2], exclude="yes")
    if(is.character(gt.as)) stop(print(gt.as))
    snps.ex <- gt.as$excluded[id %in% fsnps$id,]
    gt.as <- gt.as$keep[id %in% fsnps$id]
    if(!nrow(gt.as)) stop("No exonic snps in vcf")


    gcoord <- fread(gene.coord)
    ## Extract haps for fsnps and rsnp from reference panel ########
    if(is.numeric(snps)) {
        
        if("percentage_gc_content" %in% names(gcoord)){## newer version of gcoord from psoriasis snakefile rule geno_info using start, end, chrom, longest transcript length and GC percentage
            cis_window <- tryCatch({gcoord[gene_id==gene & chrom==chr,.(start,end)] + c(-snps, snps)},
                                   error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        } else { ## old version fp gcoord
            cis_window <- tryCatch({cl_coord(file=gene.coord,chr,gene=gene,cw=snps)}, error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        }
        
        if(is.character(cis_window)) stop(cis_window)
        if(is.list(cis_window)) cis_window <- unlist(cis_window)
        rp <- haps.range(file1=le.file,file2=h.file,cw=cis_window,population, maf)
        if(!nrow(rp)) stop("No snps extracted from the reference panel")
        rp.r=rp
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

    if(nrow(f.ase) == 0) stop("None of the fsnps are in the reference panel")
    
    rp.f=rp[rownames(rp) %in% f.ase$id, ,drop=FALSE] ## exonic snps to work with

    if(nrow(rp.f) == 0) stop("None of the fsnps are in the reference panel")
        
    ## Remove fSNPs below Fisher cut-off, if required
    if(!is.null(ex.fsnp)) {
        if(is.numeric(ex.fsnp)) {
            ## if pop!=ALL and sample.file !=NULL, select population from reference panel to check frequency of genotypes
            if(!is.null(sample.file) & population !="ALL"){
                sam <- fread(sample.file)
                ## get rows in sam  for population
                rows.panel <- sam[GROUP==population, which=T]
                ## samples in hap file are in columns and each sample corresponds to 2 cols
                ## I need to get for each sample in which col of hap file it is.
                cols <- sort(c(rows.panel*2, rows.panel*2-1))
                
                ##only use reference panel for population
                het.all <- prop_het(f.ase,rp.f[,cols, drop=F],gene)
            } else {          
                het.all <- prop_het(f.ase,rp.f,gene)
            }
            
            het.f <- het.all[pvalue > ex.fsnp,]
            if(!nrow(het.f)){
                stop("None of the fsnps are above Fisher cut-off for p-value")
            } else {
                f.ase <- f.ase[id %in% het.f$fsnp,]
                rp.f <- rp.f[rownames(rp.f) %in% f.ase$id, ,drop=FALSE]            
            }
        } else {
            f.ase <- f.ase[ !id %in% ex.fsnp,]
            rp.f <- rp.f[rownames(rp.f) %in% f.ase$id, ,drop=FALSE] 
        }
    }
    
    ## get rSNPs
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
   
    ## get ASE counts per individual for fSNPs in reference panel
    
    fsnps <- fsnps[id %in% f.ase$id, ]

    c.ase <- tot.ase_counts(x=f.ase)

    if(is.character(c.ase)) stop(c.ase)

    c.ase = filt.fsnp(c.ase,ase=1,min.ase.snp=0,n=min.ase.n, rem=NULL)
    if(is.character(c.ase)) stop(c.ase)

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
    
    
    if(nrow(fsnps))  { ### option for NB, go on for ASE

        ## select SNPs in reference panel and order c.ase as in f.ase
        c.ase=c.ase[,unlist(lapply(f.ase$id, function(i) grep(i, colnames(c.ase), value = TRUE)))]

        
        if(!is.character(c.ase)) { ## try ASE, second check
           
            ## look at ufsnps, the ones to use for ASE counts
            if(!is.null(u.esnps)){
                c.ase <- aux.in2(gene, u.esnps, c.ase, ai)
            }
           
            if(!is.character(c.ase)) { ## go on with ASE, third check   
                
                ## Check if there are sufficient hets with enough ASE counts
                c.ase = filt.fsnp(c.ase,ase=min.ase,min.ase.snp,n=min.ase.n, rem=NULL)
                
                
                if(!is.character(c.ase)) { ## go for NB only, fourth check
     
##############################  Stan inputs ###############################
                    
                    ## fsnps, pre-compute what I need because applies to every snp

                    print("Preparing stan inputs")
                    ##ai <- if(is.null(AI_estimate)) NULL else ai[id %in% f.ase$id,.(id,  NREF_post, NALT_post, Total_post, AI_post)]
                    if(!is.null(ai)) {
                        ## make all inputs with the same fSNPS
                        ai <- ai[id %in% f.ase$id,.(id, NREF_post, NALT_post,Total_post, AI_post)]
                        c.ase <- c.ase[,unlist(lapply(ai$id, function(i) grep(i, colnames(c.ase), value=T)))]
                        rp.f <- rp.f[ai$id,,drop=F]
                        f.ase <- f.ase[id %in% ai$id,]
                    }

                    print(paste("Effective number of exonic SNPs:", nrow(f.ase)))
                   

                    inp <- aux.in3(gene, ai, case=c.ase, rp.f, rp.r, f.ase, counts.g, covariates, min.ase, min.ase.n, info=info, snps.ex, prefix, out)

                    if(is.character(inp)) return(inp)

                    if(!is.null(AI_estimate)) inp <- c(inp,  min_AI=min(ai$AI_post))

                    if(!is.null(ex.fsnp) & is.numeric(ex.fsnp)) { ## add fisher test below cut-off and fisher test for all fSNPs to save
                        inp[['het.f']] <- het.f
                        inp[['het.all']] <- het.all
                    }
                    

                    
                    
                    return(list(inp=inp,model="full", probs=probs, r.tag=r.tag))
                        
                        
                    } ## NB

                }
            }
           
    }

    ## try NB

    ## use c.ase to prepare inputs with already coded functions as I will only use NB. To run NB I need fSNPs so first call to c.ase should work
    ## call again for c.ase as it may has been modified depending where code breaks for ASE

    
    c.ase <- tot.ase_counts(x=f.ase)
    inp <- aux.in3(gene, ai=NULL, case=c.ase, rp.f, rp.r, f.ase, counts.g, covariates, min.ase=0, min.ase.n=0, info=info, snps.ex, prefix, out)

    if(is.character(inp)) return(inp)
    
    return(list(inp=inp, model="NB", probs=probs, r.tag=r.tag))
    
    
} 
