## Functions to run stan scripts more efficiently
library(data.table)
library(MASS)
library(emdbook) #simulate beta binomial
##library('Matrix', lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.4");
##library('iterpc', lib.loc="/home/ev250/R/x86_64-pc-linux-gnu-library/3.3");
library(Matrix)
library(mvtnorm)
library(gridExtra)
library(ggplot2)
library(tidyr)


source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')
source('/home/ev250/Cincinatti/Functions/various.R')


#' prepare input for stan neg.only.bj.stan
#'
#' This function allows you to prepare inputs for neg.only.bj.stan
#' @param DT1 data table with total counts, 1 row for gene, only cols for counts, as in top.rsnp.22.r
#' @param DT2 data table with rsnp GT (row) for samples (colums), extra cols allowed
#' @param covar data table with covariate information, exclude intercept
#' @keywords stan input neg only
#' @export
#' @return list to input to neg.only.bj.stan
#' input.neg.only.bj()

input.neg.only.bj <- function(DT1,DT2,covar=1){

    N <- ncol(DT1)
    cov <- stan.cov(N,covar)
    y=unname(unlist(DT1))
    g.sub <- unname(unlist(DT2[,which(names(DT2) %in% paste0(names(DT1), "_GT")), with=F]))
    L <- list(N=N, K=ncol(cov)-1, Y=y, g=g.sub, cov=cov)
    return(L)
}


#' prepare input for stan negbeta.eff.stan
#'
#' This function allows you to prepare inputs for negbeta.eff.stan
#' @param x matrix with total counts, n counts, m counts, genotype and model matrix
#' @param ASE minum value for ASE counts (m counts) to allow ASE estimaltion
#' @keywords stan input trecase
#' @export
#' @return list to input to negbeta.eff.stan
#' stan.input.negbeta.eff()

stan.input.negbeta.eff <- function(x,ASE=5){
    arr.int <- x[,c(1,4)] # total counts and genotype for all samples as array of integers
   
    y <- x[x[,3]>=ASE, 2:4] # ASE counts and genotype
    v <- 0:max(y[,2]) # vector to use for sums in beta binomial likelihood stan
    K=ncol(x)-4
    mod.mat <- x[,c(5,5:ncol(x))] # model matrix as matrix
        #B <- ncol(mod.mat) # number of coefficients for neg binon regression
    #}
    
    L <- list(N=nrow(x), K=K, A=nrow(y), V=length(v), yg=arr.int, cov=mod.mat,nmg=y, v=v)
    return(L)
}


#' function for formatting input for stan
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual
#' @param g matrix of genotypes
#' @param p.hap.pairs p(H|G): output from p.h.pair
#' @param h1 matrix with haplotype 1 for n individuals
#' @param h2 matrix with haplotype 2 for n individuals
#' @param geno.exp output from tot.ase
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5, here null because I have use this function for simulations and dont without cut-off
#' @keywords stan input
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.neg.beta.prob.eff()

stan.neg.beta.prob.eff <- function(g, p.hap.pairs,h1,h2,geno.exp,ase=5, n=NULL){
   
    n.col <- grep("\\.n",names(geno.exp), value=T)
    m.col <- grep("\\.m",names(geno.exp), value=T)
    m.counts <-  rowSums(geno.exp[,m.col,with=F])
    # select ASE input when total ase counts are above threshold and 
    A <- which(m.counts>=ase)
    if(is.null(n)){
        n=1
    }
    if(nrow(geno.exp[A,][abs(rsnp)==1,])<n){ ## at least n hets with sufficient ase counts
            return("Not enough individuals with ASE counts")
            } else {
                g.trimmed <- geno.exp[A,]
                g.obs <-apply(g[A,],1,paste,collapse="") 
                h1.short <- apply(h1[A,],1,paste,collapse="")
                h2.short <- apply(h2[A,],1,paste,collapse="")
                inp <- list(yg=geno.exp[,.(y,rsnp)],gm=data.table(g.ase=g.trimmed$rsnp,m=m.counts[A]),n=list(), p=list())
                for(i in 1:nrow(g.trimmed)) {
                    p.v <- c()
                    n.v <- c()
                    ## get haps from population with p(H|G)>0
                    if(!g.obs[i] %in% colnames(p.hap.pairs)){# when the observed genotype is not compatible with the ones observed in the reference panel
                        p.v <- c(p.v,0)
                    } else {
                        haps.g <- p.hap.pairs[which(p.hap.pairs[,g.obs[i]]>0), g.obs[i]]
                                        #start with hap={h1,h2}
                                        # add n.counts, genotype associated with n is defined as in {h1,h2}
                        n.v <- c(n.v,sum(g.trimmed[i,n.col,with=F]))
                                        # get p(H|G)
                        if(length(haps.g)==1){# only one haplotype, add p(H|G)=1
                            p.v <- c(p.v,1)
                        }  else {           
                            ##match {hap1,2} with names(haps.g)
                            tmp <- strsplit(names(haps.g), ",")   
                            h1.2.p <- haps.g[which(lapply(1:length(tmp), function(x) which(h1.short[i] %in% tmp[[x]] & h2.short[i] %in% tmp[[x]]))>=1)]
                            if(length(h1.2.p)==0){# this happens when the simulated true haplotype is not in reference panel
                                p.v=0
                            } else {
                                
                                pc <- unname(h1.2.p)/sum(haps.g) # conditional p, may need to be rescaled if I have no ase counts for some of the potential haplotypes
                                
                                ##for the other haps.g get n and p
                                haps.other <- haps.g[which(names(haps.g) != names(h1.2.p))]
                                ## get swaps
                                p.v.swaps=c()  ## if I have no info for some swaps this p needs to be rescaled
                                for(j in 1:length(haps.other)) {
                                    hap.x <- strsplit(names(haps.other[j]), ",")[[1]][1]
                                    ##get the differences between obs hap and pop haps to swap ase
                                    dif <- which(unlist(strsplit(h1.short[i],split=""))!=unlist(strsplit(hap.x,split="")))
                                    if(sum(dif == ncol(g))==1){ # at least swapped rSNP
                                        dif <- which(unlist(strsplit(h2.short[i],split=""))!=unlist(strsplit(hap.x,split=""))) # swap haplotype2 to use same calculation
                                    }
                                    ## get ase: only if sum(g.trimmed[i,m.col[dif], with=F])!=0, otherwise I will be swapping snps with no counts,  equivalent to no swap.

                                    if(sum(g.trimmed[i,m.col[dif], with=F])!=0){
                                    
                                        n.v <- c(n.v,sum(g.trimmed[i,c(n.col[-dif],m.col[dif]), with=F]) - sum(g.trimmed[i,n.col[dif],with=F]))
                                        ## get p
                                        p.v.swaps <- c(p.v.swaps,unname(haps.other[j])/sum(haps.g))
                                    }
                                }
                                if(length(p.v.swaps)==0){
                                    p.v <- c(p.v,1) ## pc becomes 1 as it is the only hap with ase info
                                } else {

                                    p.v <- c(p.v, c(pc,p.v.swaps)/sum(c(pc,p.v.swaps)) )
                                }
                                
                                    
                                
                            }
                        }
                    }
                    
                    ##print("i, n.v"); print(i) ; print(n.v)   
                    inp$n[[i]] <- n.v
                    inp$p[[i]] <- p.v
                }
                ## get entries with inp$p==0, sample genotype not compatible with reference panel
                w <- which(sapply(inp$p, function(i) which(i==0))==1)
                ## remove those observations from inp$gm, inp$n and inp$p
                if(length(w)!= 0){
                    inp$gm <- inp$gm[-w,]
                    inp$p <- inp$p[-w]
                    inp$n <- inp$n[-w]
                }
                ## check if enough individuals with ase counts after removing not compatible with reference panel:
                if(nrow(inp$gm[abs(g.ase)==1,])<n){
                    return("Not enough individuals with ASE counts")
                } else {
                    return(inp)
                }
                
            }
}


#' function for formatting input for stan: removes individuals with genotypes not compatible with the reference panel
#'
#' This function allows you remove individuals from ASE with genotypes not compatible with ref panel
#' @param l list output from stan.neg.beta.prob.eff
#' @keywords stan input remove
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.rem()

stan.rem <- function(l){
    p0 <- sapply(l, function(i)  which(sapply(i$p, function(j) which(j==0))==1)) # for each simulation gives the individuals with p=0
    w <- which(sapply(p0,length)!=0) # get indexes for simulations with individuals to remove
    for(i in w){
        l[[i]]$p <- l[[i]]$p[-p0[[i]]]
        l[[i]]$n <- l[[i]]$n[-p0[[i]]]
        l[[i]]$gm <- l[[i]]$gm[-p0[[i]], ]
    }
    return(l)
}



#' function for formatting input for trecase stan, wrapping p.hap.pair and stan.neg.beta.prob.eff to avoid creating big objects when returning from p.hap.pair
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual. Starts from refernce panel input, calculates p(H|G) and then prepares stan input for full model
#' @param rp.1r numeric of haplotypes from reference panel for 1 rsnp
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param rs.hap data.table with haplotypes for the rsnp in all samples, 1 row 
#' @param geno.exp output from tot.ase for the rsnp id
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.het minimun number of het individuals for rsnp with counts >=ase, trease default 5
#' @keywords stan input
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.trecase.eff()

stan.trecase.eff <- function(rp.1r, rp.f, f.ase, rs.hap, geno.exp,min.ase=5, min.ase.het=5){
  
  rp2 <- t(rbind(rp.f,rp.1r)) ## get ref panel fsnps and 1 rsnp in the same format I have functions from simulations
  
  ## calculate P(H|G) ref panel for each rsnp                   
  tmp <- p.hap.pair.r(rp2) 
  
  ## get haps for fsnps and rsnp
  h.samp <- hap_sam(x=f.ase,y=rs.hap$id,z=rs.hap)
  g= Reduce("+", h.samp)
  
  tmp <- stan.neg.beta.prob.eff(g , p.hap.pairs=tmp, h1=h.samp[[1]], h2=h.samp[[2]], geno.exp=geno.exp, ase=min.ase, n=min.ase.het )
  
  return(tmp)
  
  
}


#' function for formatting input for trecase stan, wrapping p.hap.pair and stan.neg.beta.prob.eff to avoid creating big objects when returning from p.hap.pair,k version 2
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual. Starts from refernce panel input, calculates p(H|G) and then prepares stan input for full model
#' @param counts.g named vector  with total gene counts, names are sampleID
#' @param rp.1r numeric of haplotypes from reference panel for 1 rsnp
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param rs.hap data.table with haplotypes for the rsnp in all samples, 1 row
#' @param rec.rsnp data table with GT for rsnp in scale 0,1,-1,2 in all samples, 1 row
#' @param stan.f output from stan.fsnp.noGT.eff list with fSNPs summary data
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.het minimun number of het individuals for rsnp with counts >=ase, trease default 5
#' @keywords stan input
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.trecase.eff2()

stan.trecase.eff2 <- function(counts.g, rp.1r, rp.f, f.ase, rs.hap, rec.rsnp, stan.f,min.ase=5, min.ase.het=5) {

    ######### ASE side
    ## Start by checking if enough ind with ASE counts, otherwise abort
   
    m=stan.f$m
    ## remove "_GT" from samples in "m" in  rec.rsnp to ease computations

    setnames(rec.rsnp, sapply(names(m), function(i) grep(i,names(rec.rsnp))), names(m))
    gm <- data.table(g.ase= unname(unlist(rec.rsnp[1, names(m), with=F])) , m=m)

    ## check if enough hets

    if(sum(abs(gm$g.ase)==1) <  min.ase.het) {
        return("Not enough het individuals with sufficient ASE counts")
    }
   
  
    rp2 <- t(rbind(rp.f,rp.1r)) ## get ref panel fsnps and 1 rsnp in the same format I have functions from simulations
    
    ## calculate P(H|G) ref panel for each rsnp                   
    M <- p.hap.pair.s(rp2) 

    ## make M.cond to extract p(H|G)
    M.cond <- mat.cond(M)
    DT.cond <- data.table(summary(M.cond))  ## to easily get non-0 entries
   
    
    ## make named vector with names rownames(M) and value swapped rownames(M)
    row.split <- unlist(strsplit(rownames(M), ","))
    row.2 <- row.split[1:length(row.split) %%2 == 0]  ## take even entries
    row.1 <- row.split[1:length(row.split) %%2 ==1] ## odd entries
    s.M <- paste(row.2, row.1, sep=",")  ## invert haplotypes
    names(s.M) <-  rownames(M)
    
    ## get n counts and p for each  sample
    n.all <- list()
    pH.all <- list()
    if(any(names(stan.f) == "ai")) {
        ai.all <- list()
        vai0 <- list()
    }
    
    
    for(i in 1:length(stan.f$n)) {
        n.mat <- stan.f$n[[i]]
        ## get hap pairs for fsnps 
        hap.f.p <- lapply(colnames(n.mat), function(i) unlist(strsplit(i,",")))

        ## get rsnp GT for samples in n.mat
        gt.r <- unlist(rec.rsnp[, rownames(n.mat), with=F])
        
        ## get unique gt.n

        u.gt.r <- unique(unlist(gt.r))

        ## by gt
        u.hom <- u.gt.r[u.gt.r==0 | u.gt.r==2]
        u.het <- u.gt.r[abs(u.gt.r)==1] ## any het
        if(length(u.hom)) {

            hap.r.hom <- lapply(hap.f.p, function(i) { tmp <- sapply(u.hom, function(j) paste(paste0(i,rep(j/2,2)), collapse=","))
                names(tmp) <- paste0("g",u.hom)
                return(tmp)
            })
            
            if(length(hap.f.p)==1){
                hap.hom <- as.data.table(matrix(unlist(hap.r.hom), nrow=1, dimnames=list("",names(hap.r.hom[[1]]))))
               
            } else {
                
                hap.hom <- as.data.table(Reduce(rbind,hap.r.hom))
            }
             
        }

        if(length(u.het)){ ## first half of list for gase==1, second half for gase==-1

            hap.het <- data.table(Reduce(cbind,lapply(list(0:1, 1:0), function(i) Reduce(rbind,lapply(hap.f.p, function(j) paste(paste0(j,i, collapse=",")))))))

            names(hap.het) <- paste0("g1.", c("p","m"))
        }
 
        ## get p(H|G)
        if(length(u.hom) & length(u.het)) { ## merge first
            DT <- cbind(hap.hom,hap.het)
           
        } else {
            if(length(u.hom)){
                DT <- hap.hom
            } else {
                DT <- hap.het
            }
        }
        DT[, fhap:=colnames(n.mat)]
        keep <- grep("^g",names(DT), value=T)
        DT[, paste0("row.", keep):=lapply(keep, function(i) hap.p.no.gt2(s.M, M=M, get(i)))]
        ## assign p=0 and replace when row.g !=0
        DT[, paste0("p.", keep):=0]
        ## to add allelic imbalance
        if(any(names(stan.f) == "ai")){
            ai.mat <- stan.f$ai[[i]]
            vai.mat <-  stan.f$vai[[i]]
        } else {
            ai.mat <- NULL
            vai.mat <- NULL
        }       
        for(k in keep){

            ## sort to make sure DT.cond and row.k are in the same order
            setkeyv(DT,paste0("row.",k))
            DT[get(paste0("row.",k)) !=0, paste0("p.",k):= DT.cond[i %in% get(paste0("row.",k)), x]]
        }

        ## prepare lists for this set of individuals by GT rsnp
        set.list <- lapply(unique(gt.r), function(i) list.help(i, gt.r, n.mat, DT, m, ai.mat, vai.mat))

        ##select n and p and append n and p into 1 list respectively
        n.all <- c(n.all, Reduce("c", lapply(set.list, `[[`, 1))) 
        pH.all <- c(pH.all, Reduce("c", lapply(set.list, `[[`, 2)))
        if(any(names(stan.f) == "ai")) {
            ai.all <- c(ai.all, Reduce("c", lapply(set.list, `[[`, 3)))
            vai0 <- c(vai0, Reduce("c", lapply(set.list, `[[`, 4)))
        }

    } 
    ## n.all and pH.all need to be sorted as m

    n.all <- n.all[names(m)]
    pH.all <- pH.all[names(m)]
    if(any(names(stan.f) == "ai")) {
        ai.all <- ai.all[names(m)]
        vai0 <- vai0[names(m)]
    }
    

    ## get entries with pH.all==0, samples with  genotype of rsnp is not compatible with reference panel
    a <- sapply(pH.all, function(i) sum(i)==0)
    ## remove those observations from gm, n.all and pH.all (ai.all if applicable)
    if (any(a)){
        gm <- gm[!a,]
        n.all <- n.all[!a]
        pH.all <- pH.all[-a]
        if(any(names(stan.f) == "ai")) {
            ai.all <- ai.all[!a]
            vai0 <- vai0[!a]
        }
        
        if(nrow(gm[abs(g.ase)==1,])<min.ase.het){
            return("Not enough individuals with ASE counts")
        }
    }
    
    ## NB side
    yg <- data.table(y=unname(counts.g), rsnp=unname(unlist(rec.rsnp[, sapply(names(counts.g), function(i) grep(i,names(rec.rsnp))), with=F])))

    
    ## prepare list to return
    if(any(names(stan.f) == "ai")) {
        l <- list(yg=yg, gm=gm, n=n.all,p=pH.all, ai=ai.all, vai=vai0)
    } else {
        l <- list(yg=yg, gm=gm, n=n.all,p=pH.all)
    }
    
    return(l)
  
  
}

#' helper function to prepare list of inputs for stan.trecase.eff2
#'
#' This function allows you to format selected n counts and  pH|G based on GT
#' @param g genotype rsnp in scale 0,1,-1,2
#' @param gt.r named numeric vector genotype for rsnp in scale 0,1,-1,2 for all samples, names are sample ID
#' @param n.mat  matrix of n counts for each hap pair of fsnps concordant with reference panel (cols) and rows sample id
#' @param DT data table with hap pairs of fsnp and rsnp by genotype rsnp, and their respective probabilibities, as made in stan.trecase.eff2
#' @param m character vector with total ASE counts, names sample ID
#' @param ai.mat  matrix of allelic imbalance estimates as logit for each hap pair of fsnps concordant with reference panel (cols) and rows sample id, defaults to NULL
#' @param vai.mat  matrix of variance for allelic imbalance estimates for each hap pair of fsnps concordant with reference panel (cols) and rows sample id, defaults to NULL
#' @keywords stan GT input rsnp
#' @export
#' @return list of list of 1) list of n counts, named by sample id, 2) list of p(H|G) named by sample id 3) optional allelic imbalance named by sample id
#'
#' list.help()

list.help <- function(g,gt.r,n.mat,DT,m,ai.mat=NULL,vai.mat=NULL){
    
    sam.id <- names(gt.r)[gt.r==g]
    if(g==0 | g==2){
        n.sam <- lapply(sam.id, function(i) unname(n.mat[i,]))     
        pH.sam <- lapply(sam.id, function(i) unname(unlist(DT[, paste0("p.g",g), with=F])))
        if(!is.null(ai.mat)){
            ai.sam <- lapply(sam.id, function(i) unname(ai.mat[i,]))
            vai.sam <- lapply(sam.id, function(i) unname(vai.mat[i,]))
        }        
    }
    if(g==1 | g==-1){ ## swap haps of fsnps to avoid changing g.ase
        mn.mat <- apply(n.mat[sam.id,,drop=F], 2, function(i) m[sam.id]-i)
        if(!is.matrix(mn.mat)){
            mn.mat <- matrix(mn.mat, nrow=1, dimnames=list("",names(mn.mat)))
        }          
        mat2 <- cbind(n.mat[sam.id,,drop=F],mn.mat)
        n.sam <- lapply(1:nrow(mat2), function(i) unname(mat2[i,]))
        if(!is.null(ai.mat)){
            ## swap ai (logit to -logit)
            ai.swap.mat <- apply(ai.mat[sam.id,,drop=F], 2, function(i) -i)
            if(!is.matrix(ai.swap.mat)){
                ai.swap.mat <- matrix(ai.swap.mat, nrow=1, dimnames=list("",names(ai.swap.mat)))
            }          
            matAI <- cbind(ai.mat[sam.id,,drop=F],ai.swap.mat)
            ai.sam <- lapply(1:nrow(matAI), function(i) unname(matAI[i,]))
            matVAI <- cbind(vai.mat[sam.id,,drop=F], vai.mat[sam.id,,drop=F])## no need to swap
            vai.sam <- lapply(1:nrow(matVAI), function(i) unname(matVAI[i,])) 
        }
        
        p <- c("p.g1.p", "p.g1.m")
        if(g==-1){
            p <-p[length(p):1]
        }
        ## order DT$fhap as names(n.mat)
        DT <- DT[order(match(fhap, colnames(n.mat)))]
        pH.sam <- lapply(sam.id, function(i) unname(unlist(DT[, p, with=F])))
    }

    names(n.sam) <- names(pH.sam) <- sam.id
    if(!is.null(ai.mat)){
        names(ai.sam) <- names(vai.sam) <- sam.id
        l <- list(n.sam, pH.sam, ai.sam, vai.sam)
    } else {
        l <- list(n.sam, pH.sam)
    }

    return(l)
}




#' sub-function to deal with covariates when preparing input for stan 
#'
#' This function allows you to deal with covariates for stan input
#' @param N, number of individuals with covariate information
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input covariates
#' @export
#' @return matrix with covariates compatible with stan: two cols with ones plus covs, first col will be ignored by stan but allows to enter data as matrix avoidind vector when no covs
#' stan.cov()

stan.cov <- function(N, covar=1){

    cov <- matrix(1,nrow=N, ncol=2) # I need at least 2 cols so stan will recognise covar as a matrix, which helps with downstream calculations in stan. first extra column will be ignored
    if((is.numeric(covar) & length(covar)==N) | is.matrix(covar)){
        cov <- cbind(cov, covar) ## add covariates   
    }
    rownames(cov) <- NULL
    return(cov)
}



#' prepare input for stan neg.beta.prob.priors.eff.phasing.stan
#'
#' This function allows you to prepare inputs for negbinombeta.ase.prob.phasing.stan
#' @param x element from output list from stan.neg.beta.prob.eff
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase haplotype probabilities 
#' @export
#' @return list to input to stan neg.beta.ase.prob.priors.eff.phasing.stan
#' in.neg.beta.prob.eff()

in.neg.beta.prob.eff <- function(x,covar=1){
    
    N <- nrow(x[[1]]) # number of individuals
    cov <- stan.cov(N,covar)
    # n. counts    
    n.v <- do.call(c,unname(x$n))
    s <- unname(sapply(x$n,length))
    L <- sum(s)
    #p(H)
    p.v <- do.call(c,x$p)
    # v
    v <- 0: max(x$gm$m)
    LL=list(N=N, A=nrow(x$gm), L=L, K=ncol(cov)-1, M=length(v), Y=x$yg$y, g=x$yg$rsnp, gase= x$gm$g.ase, m=x$gm$m, n=n.v, pH=p.v, s=s, v=v, cov=cov)
    return(LL)
}


#' prepare input for stan neg.beta.prob.priors.eff2.stan
#'
#' This function allows you to prepare inputs for neg.beta.prob.phasing.priors.eff.test.stan
#' @param x element from output list from stan.neg.beta.prob.eff
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase haplotype probabilities 
#' @export
#' @return list to input to stan neg.beta.ase.prob.priors.eff2.stan
#' in.neg.beta.prob.eff2()

in.neg.beta.prob.eff2 <- function(x,covar=1){
    
    N <- nrow(x[[1]]) # number of individuals
    cov <- stan.cov(N,covar)
    ## n. counts    
    n.v <- do.call(c,unname(x$n))
    s <- unname(sapply(x$n,length))
    L <- sum(s)
    ## p(H)
    p.v <- do.call(c,unname(x$p))

    if(any(names(x) == "ai")){
        ai.v <- do.call(c, unname(x$ai))
        vai.v <- do.call(c, unname(x$vai))
        LL=list(N=N, A=nrow(x$gm), L=L, K=ncol(cov)-1, Y=x$yg$y, g=x$yg$rsnp, gase= x$gm$g.ase, m=x$gm$m, n=n.v, pH=p.v, ai0=ai.v, sdai0=sqrt(vai.v) , s=s, cov=cov)
    } else {
        LL=list(N=N, A=nrow(x$gm), L=L, K=ncol(cov)-1, Y=x$yg$y, g=x$yg$rsnp, gase= x$gm$g.ase, m=x$gm$m, n=n.v, pH=p.v, s=s, cov=cov)
    }
    return(LL)
}





#' Sub function for simulating haplotypes with ASE, fSNPs plus one rSNP: total and ase expression allowing subgroups within ASE
#'
#' This function allows you to simulate total and ASE for haplotype simulation
#' @param f number of phasing snps
#' @param g vector of genotypes
#' @param betas vector with intercept and coefficients for negative binomial model matrix
#' @param DT1 data table  with intercpept (first col 1's) and covariates for calculating the mean for negative binomial
#' @param phi negative binomial overdispersion
#' @param theta beta binomial overdispersion
#' @param p prob for ASE in het individuals, connected to bj, p=exp(bj)/(1+exp(bj))
#' @param f.ase ASE fraction relative to total counts
#' @param betas2 regression coefficients for subgroup analysis of genotype effect on counts
#' @param DT2  data table with intercpept (first col 1's) and covariates for calculating the mean for ASE
#' @keywords simulation negative binomial trec counts and beta binomial ASE counts
#' @export
#' @return matrix with counts, hap.a counts, total ASE counts  and genotype
#' tot.ase.sub()

tot.ase.sub <- function(f,g,betas,DT1, phi,theta,p,f.ase,betas2,DT2){
    bj=log(p/(1-p))
    lmu=as.matrix(DT1) %*% betas
    DT1[,rsnp:=g[,ncol(g)]]
    DT1[,y:=rnbinom(1,size=phi,mu=exp(lmu))] # default rsnp=0
    DT1[rsnp==1, y:=rnbinom(1,size=phi, mu=exp(lmu + log(1+ exp(bj)) -log(2)))]
    DT1[rsnp==2, y:=rnbinom(sum(rsnp==2),size=phi, mu=exp(lmu + bj))]
    
    ##############################################################################
### working here, adapting tot.ase to incorporate extra features in teh function #####################################################################################
    
    # calculate ASE for fSNPs=1 based on rSNP
    DT <- cbind(DT,g[,1:f])
    names(DT)[3:(f+2)] <- paste0("fsnp.",1:f)
    for(i in  paste0("fsnp.",1:f)){
        DT[, paste0(i,".n"):=vector(mode="integer",nrow(DT))][,paste0(i,".m"):=vector(mode="integer",nrow(DT))]
        
        DT[rsnp!=1 & get(i)==1, paste0(i,".n"):=rbetabinom(length(which(DT$rsnp!=1 & DT[,i,with=F]==1)), prob=0.5, size=round(f.ase*y), theta)]
        
        DT[rsnp==1 & get(i)==1, paste0(i,".n"):=rbetabinom(length(which(DT$rsnp==1 & DT[,i,with=F]==1)), prob=p, size=round(f.ase*y), theta)]
        
        DT[get(i)==1, paste0(i,".m"):=as.integer(round(f.ase*y))]
        }
        return(DT)
}




#' fix haplotypes and prepare input for stan neg.beta.prob.priors.eff.phasing.stan
#'
#' This function allows you to prepare fix haplotypes (as trecase) and prepare input for negbinombeta.ase.prob.phasing.stan
#' @param x list output from stan.neg.beta.prob.eff or from sim.haps with option out=in.neg.beta.prob.eff
#' @keywords stan input fixed haplotype
#' @export
#' @return list to input to stan neg.beta.ase.prob.priors.eff.phasing.stan
#' fixhap.eff()

fixhap.eff <- function(x){
    #n.fix1 <- lapply(seq_along(x), function (i) lapply(x[[i]]$n, `[[`, 1)) # extract first element in "n"
    x$n <- lapply(x$n, `[[`, 1) # repalce x$n with the first element in "n"
    x$p <- rep(list(1), length(x$n)) # replace x$p with list of 1's

   fix.input1 <- in.neg.beta.prob.eff(x=x,covar=1)
    return(fix.input1)
}

####################### help with setting simulation conditions ######################


#' simulate population, and from the population take a sample. Calculate the difference between p(H) pop and p(H) sample
#'
#' This function allows you to choose simulation parameters according to the difference in p(H) pop and p(H) sample
#' @param f number of fSNPs
#' @param m mean of each SNP (f+1) to simulate haplotypes using multinormal distribution
#' @param cov covariance between all SNPs (fSNPs and 1 rSNP)
#' @param N number of individuals in population
#' @param n number of individuals to sample
#' @param r2 whether to return r2
#' @param maf whether to return maf
#' @keywords simulations p(H) 
#' @export
#' @return sum of the square differences between p(H) pop and sample
#' pH.diff()

ph.diff <- function(f=3,m,cov,N,n,r2=NULL,maf=NULL){

    ## pop
    hm1 <- sim.pop.haps(f=3,m,cov,N,r2=TRUE,maf=TRUE) # N number of individuals in population
    ##samples
    s1 <- hm1[[1]][sample(1:nrow(hm1[[1]]),size=n,replace=TRUE),]
    s2 <- hm1[[1]][sample(1:nrow(hm1[[1]]),size=n,replace=TRUE),]
    ## calculate P(H) pop
    pop.haps <- data.table(apply(hm1[[1]],1,paste,collapse=""))
    pop.freq <- pop.haps[,.(N=.N/nrow(pop.haps)),by=pop.haps]
    ##pop.freq <- data.table(prop.table(table(pop.haps)))
    ## calculate P(H) samples
    s.haps <- data.table(apply(rbind(s1,s2), 1, paste,collapse=""))
    samp.freq <- s.haps[,.(N = .N/nrow(s.haps)), by=s.haps]
    tmp <- merge(pop.freq,samp.freq,by="V1", all=T)
    tmp[is.na(N.y), N.y:=0]
    tmp[,dif.sq:= (N.x- N.y)^2 /N.x]
    tmp[,dif:=abs(N.x- N.y) /N.x]
    return(sum(tmp$dif.sq))
}

#'  Calculate the square difference between p(H) pop and p(H) sample
#'
#' This function allows you calculate the square difference between p(H) pop and p(H) sample starting from the observed haplotypes
#' @param x matrix of haplotypes for population
#' @param y matrix of haplotypes for sample
#' @keywords haplotype frequencies differences
#' @export
#' @return sum of the square differences between p(H) pop and sample
#' pH.diff2()

ph.diff2 <- function(x,y){
     ## calculate P(H) pop
    pop.haps <- data.table(apply(x,1,paste,collapse=""))
    diff <- ph.diff3(pop.haps,y)
    return(diff)
}
                              


#'  Calculate the square difference between p(H) pop and p(H) sample
#'
#' This function allows you calculate the square difference between p(H) pop and p(H) sample starting from the observed haplotypes and the population already formatted
#' @param x data table of named haplotypes for population (rows = #individuals *2, 1 col)
#' @param y matrix of haplotypes for sample
#' @param z cut-off for population hap frequency, only compute the sqd for haplotypes which pop frequency is above cut-off
#' @keywords haplotype frequencies differences
#' @export
#' @return sum of the square differences between p(H) pop and sample
#' pH.diff3()

ph.diff3 <- function(x,y, z=NULL){
    ## calculate P(H) pop
    pop.freq <- x[,.(N=.N/nrow(x)),by=V1]
    if(!is.null(z)){
        pop.freq <- pop.freq[N>z,]
    }
    ## calculate P(H) samples
    s.haps <- data.table(apply(y, 1, paste,collapse=""))
    samp.freq <- s.haps[,.(N = .N/nrow(s.haps)), by=V1]
    tmp <- merge(pop.freq,samp.freq,by="V1", all.x=T)
    tmp[is.na(N.y), N.y:=0]
    tmp[,dif.sq:= (N.x- N.y)^2 /N.x]
    #tmp[,dif:=abs(N.x- N.y)/N.x ]
    return(sum(tmp$dif.sq))
}


##################################################################################################
## Functions for no GT in rsnp, getting it from rp with p(G)


#' sub-function for formatting input for stan with no GT for rsnp
#'
#' This function allows you to get p(H) entries compatible with hap1|hap2 of fsnps and 1 rsnp 
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param v character vector with haplotype pairs (separeated by ",") of fsnps and all possible GT for rsnp
#' @keywords stan pH unknown rsnp genotype
#' @export
#' @return named vector with pH, names correspond to paste(v , GT, sep=".")
#'
#' pH.no.gt()

pH.no.gt <- function(M,v){
    p <- lapply(v, function(i) {
        tmp <- M[i,]
        pi <- M[i, M[i,]>0]
        names(pi) <-  names(tmp)[which(tmp>0)]
        return(pi)
    })
    names(p) <- v
    p <- unlist(p)
    return(p)
}


#' sub-function for formatting input for stan with no GT for rsnp
#'
#' This function allows you to check whether hap.pairs (rownames of matrix M)  are in a input vector allowing for swap 
#' @param M matrix  with hap pairs and GT prob for reference panel
#' @param y character vector with haplotype pairs (separeated by ",") to match in x (allowing swaps)
#' @keywords stan pH unknown rsnp genotype
#' @export
#' @return character vector with names y and value the rownumber in M, if no match row=0
#'
#' hap.p.no.gt()

hap.p.no.gt <- function(M,y){
    s.M <- sapply(rownames(M), function(i) paste(rev(unlist(strsplit(i, ","))),collapse=",")) ## named vector with names rownames(M) and value swapped rownames(M)
    ## for each element of y get the row number in M either swap or unswap
    l <- lapply(y, function(i) unique(c(which(names(s.M)==i) , which(s.M==i))))
    l[which(lapply(l,length) ==0)] <- 0 ## replace  no matches with 0
    names(l) <-  y
    ul <- unlist(l)
    return(ul)
}


#' sub-function for formatting input for stan with no GT for rsnp, similar to hap.p.no.gt but avoids re-doing fisrt step (s.M)
#'
#' This function allows you to check whether hap.pairs (rownames of matrix M)  are in a input vector allowing for swap
#' @param s.M named vector with names rownames(M) and value swapped rownames(M)
#' @param M matrix  with hap pairs and GT prob for reference panel
#' @param y character vector with haplotype pairs (separeated by ",") to match in x (allowing swaps)
#' @keywords stan pH unknown rsnp genotype
#' @export
#' @return character vector with names y and value the rownumber in M, if no match row=0
#'
#' hap.p.no.gt2()

hap.p.no.gt2 <- function(s.M, M,y){
    
    ## for each element of y get the row number in M either swap or unswap
    l <- lapply(y, function(i) unique(c(which(names(s.M)==i) , which(s.M==i))))
    l[which(lapply(l,length) ==0)] <- 0 ## replace  no matches with 0
    names(l) <-  y
    ul <- unlist(l)
    return(ul)
}


#' sub-function for formatting input for stan with no GT for rsnp
#'
#' This function allows you to get p(H), n and GT for rsnp for non-observed hap of fsnps and rsnp compatibles with observed GT for fsnps and haplotypes in reference panel
#' @param M matrix with p(H) row hap.pairs and cols GT for fsnps plus rsnp from reference panel
#' @param M.cond  matrix with p(H|GT fsnps), row hap.pairs and cols GT for fsnps from reference panel, output from mat.col
#' @param m.trim  matrix with ASE information
#' @param i index for selecting rows of m.trim
#' @param n.col names for cols in mat with n counts
#' @param m.col names for cols in mat with m counts
#' @param DT data table with cols: GT.ob (observed GT,h1.2, hap1,hap2 and rows.comp row in M that matches h1.2, constructed in stan.ase.no.gt
#' @param j GT.ob, to work within loop
#' @keywords stan input unknown rsnp genotype
#' @export
#' @return list with partial inputs for stan p= collection of p(H), n=n counts, r= gt for rsnp in scale 0,1,-1,2
#'
#' sw.no.gt()

sw.no.gt <-function(M,M.cond,DT,m.trim,i,j,n.col,m.col){
    ##collect terms in vectors
    p.v <- c()
    n.v <- list()
    r.v <- c()
    gt.f <- substr(j,1,nchar(j)-1)
    ## select rows to exclude from Mcond becuase they are the obs hap
    DT.sub <- DT[GT.ob==j,rows.comp]
    ex <- which(DT.sub != rep(0,length(DT.sub)))  ## position in DT.sub to exclude
    ## select swaps for GT
    if(length(ex) == 0){ ## no rows to exclude
        hap.others <- which(M[,j]>0)
    } else {
        hap.others <- which(M[-DT.sub[ex], j] > 0)
    }
    s <-sapply(hap.others, function(i) length(i) >0) 
    if(length(s)>0){
        hap.others <- names(hap.others[which(s)])
        
        gt <-  as.numeric(substr(j, nchar(j), nchar(j))) ## gt  rsnp to swap
        ## get hap1 for GT in sample, first entry if more than 1
        hap1 <-unlist(strsplit( DT[GT.ob==j,][1,h1.2], ","))[1]
        ## get haps to swap
        hap.x <- lapply(hap.others, function(i) unlist(strsplit(i,","))[1])  ## hap 1
        for(k in 1:length(hap.x)){
            ##get the differences between hapx hap1 to swap ase
            dif <- which(unlist(strsplit(hap1,split=""))!=unlist(strsplit(hap.x[[k]],split="")))
            ## remove dif entry corresponding to rsnp
            wr <- which(dif==nchar(hap.x[[k]]))
            if(length(wr)!=0){
                dif <- dif[-wr]
            }
                                        
            ## get ase: I remove the condition sum(g.trimmed[i,m.col[dif], with=F])!=0, which avoids swapping snps with no counts,because I need to accout for this hap in p.
            ##if(sum(m.trim[i,m.col[dif]])!=0){
            n.v[[k]] <-rowSums(m.trim[i,c(n.col[-dif],m.col[dif]), drop=FALSE]) - rowSums(m.trim[i,n.col[dif], drop=FALSE])
            ## get p
            p.v <- c(p.v, M.cond[hap.others[k],as.character(gt.f)]) ## conditional p
            ## get GT rsnp coded as 0,1,-1,2
            if(gt==0 | gt==2){
                gt.r <- gt
            } else {
                h1.rsp <- gt-as.numeric(substr(hap.x[k], nchar(hap.x[k]), nchar(hap.x[k]))) ## 
                gt.r <- ifelse(h1.rsp==1,1,-1)
            }
            r.v <- c(r.v,gt.r)
            ##}
        }
    }
    l <- list(p=p.v,n=Reduce(cbind,n.v),r=r.v)
    
    return(l)    
}


#' function for formatting input for stan for ASE model with no GT for rsnp
#'
#' This function allows you to format input for stan ASE model per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param l list with matrix1 hap1 for fsnps row samples and matrix2 same but hap2
#' @param m matrix with ase counts per fsnp, output from tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals het for rsnp with counts >=ase, trease default 5, here null because I have used this function for simulations without cut-off
#' @param gf.comp named vector with GT for fsnp for samples that pass ase and n cut-offs, output from sel.ind.no.gt, defaults to null as it can be computed with input data but avoids recalculation when calling this function from stan.full.no.gt
#' @param M.cond matrix with conditional prob for an haplotype pair of fsnps and rsnp given the GT of the fsnps, output from mat.col. Defaults to null as it can be computed with input data but avoids recalculation when calling this function from stan.full.no.gt
#' @keywords stan input ase unknown rsnp genotype
#' @export
#' @return list of 1)  total ase counts (m). 2) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 3) list of vectors p: each vector with p for each haplotype corresponding to n, 4) list of vectors g: each vector g for the rsnp.
#'
#' stan.ase.no.gt()

stan.ase.no.gt <- function(M,l,m,ase=5, n=NULL,gf.comp=NULL, M.cond=NULL){
    ## get GT for fsnps in samples with sufficient ASE at least in n individuals and compatible with ref panel
    if(is.null(gf.comp)){
        gf.comp <- sel.ind.no.gt(M,l,m,ase,n)
    }
    if(length(gf.comp) == 1){ ## error message from sel.ind.no.gt
        return(gf.comp)
    } 
    m.col <- grep("\\.m",colnames(m), value=T)
    n.col <- grep("\\.n",colnames(m), value=T)
    m.trim <- m[names(gf.comp),,drop=F] ## only inds with ASE info and GT compatible with ref panel
    if(nchar(gf.comp[1])==1){ ## one fsnp, tested in first element of gf.comp
        h.short <- lapply(l, function(i) i[names(gf.comp),])
        m.counts <- m.trim[,m.col]
    } else {
        h.short <- lapply(l, function(i) apply(i[names(gf.comp),,drop=F],1,paste,collapse=""))
        m.counts <-  rowSums(m.trim[,m.col, drop=F])
    }
    
    ## get matrix p(H|GT fsnps)
    if(is.null(M.cond)){
        M.cond <- mat.col(M)
    }
    inp <- list(m=m.counts,n=list(), p=list(),g=list())
    ## to make it faster: get inputs for unique(haplotypes) and then expand it to all samples
    hap.s <- Reduce(cbind,h.short)
    colnames(hap.s) <- c("hap1", "hap2")
    hap.dt <- data.table(hap.s,keep.rownames=T)
    u.haps <- unique(hap.s,MARGIN=1)
    u.haps.dt <- data.table(u.haps)
    ## add geno fsnps
    if(length(m.col)==1) { ## 1 fsnp
        u.haps.dt[,g.fsnps:=hap1+hap2]
    } else {
        u.haps.dt[,g.fsnps:=Map(add.geno,hap1,hap2)]
    }
    
    ## get indexes for u.haps in hap.dt
    if(nrow(u.haps.dt)==1){## general form doesnt work only accepts one value
        u.haps.dt[, ind:= paste(which(hap.dt$hap1==hap1 & hap.dt$hap2==hap2), collapse=",")]
        idx=as.numeric(unlist(strsplit(u.haps.dt$ind, ",")))
    } else {
            
        u.haps.dt[, ind:=lapply(1:nrow(u.haps.dt), function(i) which(hap.dt$hap1==hap1[i] & hap.dt$hap2==hap2[i]))]
    }    
    for(i in 1:nrow(u.haps.dt)) {
        p.v <- c()
        if(nrow(u.haps.dt)>1){
            idx <- unlist(u.haps.dt$ind[i]) ## index individuals with same hap combinations
        }
        
        n.v.mat <- matrix(,nrow=length(idx), ncol=0)  ## for all individuals with same haps, each row is the n for each individual, cols will be all possible haps. ncol=length(p.v)=length(r.v)
        r.v <- c()  ## to collect GT of rsnp in scale 0,1,-1,2 with 1= 0|1 and -1=1|0
        ## add rsnp, any GT to each hap
        z <- lapply(u.haps[i,], function(j) paste0(j,0:1))
        ## get all possible haps pairs for hap1 and hap2 (fsnps and rsnp)
        h1.2 <- as.vector(outer(z$hap1,z$hap2, paste, sep=","))
        ## get rows in M compatible with obs haps
        rows.comp <- hap.p.no.gt(M, y=h1.2)
        ## get GT for h1.2
        GT.ob <- as.vector(outer(z$hap1,z$hap2, add.geno)) ## to keep same order as h1.2
        ## summarize h1.2, with genotype with rows compatible with haps in M into 1 object
        sum.gt <- data.table(GT.ob,h1.2,rows.comp)
        ## get pop entries compatibles with GT
        g.col <- which(GT.ob %in% colnames(M))                      
        ## select entries in sum.gt with GT compatible with reference panel, some rSNP GT may not be compatible
        sum.gt <- sum.gt[g.col,]
        ## in each col work out obs hap if present and swaps
        ## get GT fsnps
        gt.f <- u.haps.dt$g.fsnps[i]
        for(j in unique(sum.gt$GT.ob)) {
            ## check if obs hap
            if(sum(sum.gt[GT.ob==j,rows.comp]) !=0) { ## at least one obs hap
                r <- which(sum.gt[GT.ob==j,rows.comp] !=0)  ## get the one(s) in ref panel
                pH1 <- M.cond[sum.gt[GT.ob==j,rows.comp][r] , as.character(gt.f)] ## get pH for obs hap, pop hap1|hap2 may be swapped relative to mine but p(H) is the same
                if(nrow(n.v.mat)>1){
                    n.v.mat <- cbind(n.v.mat, replicate(length(pH1), rowSums(m.trim[idx,n.col,drop=FALSE])))
                } else {
                    n.v.mat <- cbind(n.v.mat, t(replicate(length(pH1), rowSums(m.trim[idx,n.col,drop=FALSE]))))
                }
                
                ##n.v <- c(n.v, rep(sum(m.trim[i,n.col]), length(pH1))) ## repeat n for each possible g of rsnp, n is the same.
                p.v <- c(p.v,unname(pH1))
                ## get GT by hap1,2 for rsnp to code for -1 and 1
                if(length(sum.gt[GT.ob==j,h1.2][r])==1){## only 1 rsnp GT
                    gt.r <- as.numeric(unname(sapply(strsplit(sum.gt[GT.ob==j,h1.2][r], ","), function(i) substr(i, nchar(i), nchar(i)))))
                    if(sum(gt.r)!=1){
                        gt.r <- sum(gt.r)
                    } else {
                        h1.r <- sum(gt.r)-gt.r[1]
                        gt.r <- ifelse(h1.r==0, -1, 1)
                    }
                    
                } else { ## each row is 1 possible GT for rsnp
                    gt.r <- data.table(apply(sapply(strsplit(sum.gt[GT.ob==j,h1.2], ","), function(i) substr(i, nchar(i), nchar(i))), 1, as.numeric))
                    gt.r$gt <- rowSums(gt.r)
                    gt.r[, rec:=ifelse(gt-V1==0,-1,1)]
                    gt.r <- gt.r$rec
                }
                
                r.v <- unname(c(r.v,gt.r))
            }                          
            ## check if there are other haps compatible with GT of fsnps
            in.stan1 <- sw.no.gt(M, M.cond, DT=sum.gt, m.trim,i=idx,j=j,n.col,m.col)
            if(all(!sapply(in.stan1, is.null))) { 
                
                p.v <- c(p.v,in.stan1$p)
                ##n.v <- c(n.v,in.stan1$n)
                r.v <- c(r.v,in.stan1$r)
                if(nrow(n.v.mat)==1){
                    n.v.mat <- cbind(n.v.mat, in.stan1$n)
                } else {
                    if(!is.matrix(in.stan1$n)){ ## make matrix if numeric
                        
                        n.v.mat <- cbind(n.v.mat,matrix(in.stan1$n, ncol=1,dimnames=list(names(in.stan1$n))))
                    } else {
                        n.v.mat <- cbind(n.v.mat,in.stan1$n)
                    }
                    
                }
            }
            
        }


    ##print("i, n.v"); print(i) ; print(n.v)   
    inp$n[idx] <- lapply(1:nrow(n.v.mat), function(i) unname(n.v.mat[i,]))
    inp$p[idx] <- rep(list(p.v), length(idx))
    inp$g[idx] <- rep(list(r.v), length(idx))
    ##print(i)
    }
    
    return(inp)        
}


#' Function for making matrix of P(H snps|GT snps) (colSums = 1) starting from matrix p(H)
#'
#' conditional matrix by cols, divides each col of M by its sum
#' @param M matrix
#' @keywords matrix conditional
#' @export
#' @return matrix
#' mat.cond()

mat.cond <- function(M){
    
    v=1/colSums(M)
    diag <- as(Diagonal(length(v), v), "sparseMatrix")
    tmp <- M %*% diag
    return(tmp)
}
       
#' Function for making matrix of P(h|GT fSNPs) (colSums = 1) starting from matrix of h vs GT (fsnps+rsnp)  
#'
#' collapsing matrix by columns 
#' @param M matrix
#' @keywords matrix collapsing
#' @export
#' @return matrix
#' mat.col()

mat.col <- function(M){
    x <- as.data.table(summary(M))
    x[, r.names:=row.names(M)[i]]
    x[,c.names:=colnames(M)[j]]
    n <- nchar(x$c.names[1])
    x[,newcols:=as.character(substr(c.names,1, n-1))]
    x <- x[,.(value=sum(x)),by=c("r.names","newcols")]
#######  make sparse matrix, get index for newcols (j)
    ## sort by newcols so duplicated cols are together and then get number of repetitions using rle, but need to keep col.order to get back right order before return
    col.order <- unique(x$newcols)
    setkey(x,newcols,r.names)
    index.j <- rle(x$newcols)
    x[, j:=rep(1:length(index.j$value), index.j$length)]
    ## need to reorder x as in M, r.names same order as rownames(M), and colnames as col.order
    tmp <- sparseMatrix(i=1:nrow(x), j=x$j, x=x$value, dimnames=list(x$r.names, index.j$value))
    tmp <- tmp[match(rownames(M),rownames(tmp)), match(col.order,colnames(tmp))]
    ##v=1/colSums(tmp)
    v=1/apply(tmp,2,sum)
    diag <- as(Diagonal(length(v), v), "sparseMatrix")
    tmp <- tmp %*% diag
    return(tmp)
}

## old version, slower
mat.col2 <- function(M){
    GTf <- sort(sapply(colnames(M), function(i) substr(i, 1, (nchar(i)-1) )))
    u=unique(GTf)
    tmp <- Matrix(0,nrow=nrow(M),ncol=length(u), dimnames=list(rownames(M),u))
    for(i in u){
        n <- names(GTf)[GTf==i]
        if(length(n)==1){
            tmp[,i] <- M[,n]
            } else {
                tmp[,i] <- rowSums(M[,n])
            }
    }
    v=1/colSums(tmp)
    tmp <- mat.scalar(v,tmp)
    return(tmp)
}


#' function for formatting input for stan for FULL  model with no GT for rsnp
#'
#' This function allows you to format input for stan FULL MODEL per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual
#' @param counts data table  with gene counts per sample for a gene
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param l list with matrix1 hap1 for fsnps (rownames samples, cols GT for each fsnp) and matrix2 same but hap2
#' @param m matrix with ase counts per fsnp, output from tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals  with counts >=ase, trease default 5.
#' @param rna whether to allow for missing values in GT, used when genotyping by rna, defaults to NULL
#' @keywords stan input ase unknown rsnp genotype
#' @export
#' @return list of 2 elements for samples with sufficient ase counts: 1) NB, list  with first element total counts for gene and second element list of p(GT of rsnp|GT fsnp)>0 for GT rsnp 0,1,2. 2) ase: list with ase counts total ase counts (m). 2) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 3) list of vectors p: each vector with p for each haplotype corresponding to n, 4) list of vectors g: each vector g for the rsnp.
#'
#' stan.full.no.gt()

stan.full.no.gt <- function(counts, M,l,m,ase=5, n=5,rna=NULL){
    ## get GT for fsnps in samples with sufficient ASE at least in n individuals and compatible with ref panel
    gf.comp <- sel.ind.no.gt(M,l,m,ase,n)
    if(is.null(rna)){
        if(is.null(names(gf.comp))) { ## error message from sel.ind.no.gt
            return(gf.comp)
        }   
    }    
    ## get input for NB, consider all individuals irrespective of ASE counts
    ## get matrix p(H|GT fsnps)
    M.cond <- mat.col(M)
    ## make vector with geno for the corresponding hap.pair from rownames(M) for reference panel
    gt.rf <- Reduce(add.geno, lapply(1:2, function(i) sapply(strsplit(rownames(M), ","), `[[`, i)))
    ## get genotypes of fsnps for all samples
    gf.all <- sel.ind.no.gt(M,l,m,ase=0,n=0)
    if(all(gf.all =="GT fsnps is not compatible with Ref panel" )) { ## error message from sel.ind.no.gt for all samples
        return(gf.all)
    }
    if(any(gf.all =="GT fsnps is not compatible with Ref panel" )) {## error message from sel.ind.no.gt for all samples
        gf.all <- gf.all[gf.all !="GT fsnps is not compatible with Ref panel"]
    }
    
    ## for each gf.all get p(gf+rsnp | gf)
    u.gf.all <- as.character(unique(gf.all))
    p.g <- lapply(u.gf.all, function(i) {
        w <-  which(M.cond[,i]!=0)
        tmp <- M.cond[w,i]
        names(tmp) <- substr(gt.rf[w],nchar(gt.rf[1]), nchar(gt.rf[1])) ## last character from gt.rf is geno rsnp
        tmp <- tapply(tmp, names(tmp), sum) ## add elements in vector with same genotype
        return(tmp)
    })
    names(p.g) <- u.gf.all
    ## assign p.g to each sample in gf.all
    p.g.comp <- rep(list(p.g[[1]]), length(gf.all)) ## start with first element and replace
    ##replace if necessary:
    if(length(u.gf.all) >1){
        for(i in 2:length(u.gf.all)){
            w <- unname(which(gf.all == u.gf.all[i]))
            p.g.comp[w] <- rep(list(p.g[[i]]), length(w))
        }
    }
    
    ##  total counts for relevant samples
    y <- counts[,names(gf.all),with=F]
    
    ## get input for ASE
    if(is.null(names(gf.comp))){ ## for rna only, to allow  NB side to be returned
        ase.in="No ASE for these samples, likely homo fsnps"
    } else {
        ase.in <- stan.ase.no.gt(M,l,m,ase,n,gf.comp,M.cond)
    }
    return(list(NB=list(counts=y,p.g=p.g.comp),ase=ase.in))              
}



#' sub-function for formatting input for stan  with no GT for rsnp
#'
#' This function allows you to select individuals with sufficent ASE information for FULL or ASE model with rsnp no GT
#' @param M  matrix with p(H), row hap.pairs and cols GT for fsnps and rsnp from reference panel, output from p.hap.pair
#' @param l list with matrix1 hap1 for fsnps row samples and matrix2 same but hap2
#' @param m matrix with ase counts per fsnp, output from tot.ase_counts
#' @param ase ase cut-off default is 5 counts, trecase default 5
#' @param n minimun number of individuals with sufficient ASE counts trease default 5
#' @keywords stan input ase unknown rsnp genotype
#' @export
#' @return named vector with gt for fsnps, names are individual ID
#'
#' sel.ind.no.gt()

sel.ind.no.gt <- function(M,l,m,ase,n){
    m.col <- grep("\\.m",colnames(m), value=T)
    m.counts <-  rowSums(m[,m.col, drop=FALSE])
    
    ## select ASE input when total ase counts are above threshold
    A <- which(m.counts>=ase)
    if(length(A)<n){ ## at least n individuals with sufficient ase counts (unknown genotype)
        return("Not enough individuals with ASE counts")
    } 
    g <- Reduce("+",l)
    g.obs <-apply(g[A,,drop=F],1,paste,collapse="")    
    names(g.obs)=names(A)
    ## get GT for fsnps in ref panel
    gt.fsnps <- sapply(colnames(M), function(i) substr(i,1,nchar(i)-1))
    ## select inds that have GT fsnps compatible with ref panel
    g.comp <- g.obs[g.obs %in% gt.fsnps]
    if(length(g.comp) == 0 ){
        return("GT fsnps is not compatible with Ref panel")
    } else {
        return(g.comp)
    }
    
}



#' prepare input for stan neg.beta.noGT.rsnp.priors.eff.stan
#'
#' This function allows you to prepare inputs for stan neg.beta.noGT.rsnp.priors.eff.stan
#' @param x output list from stan.full.no.gt
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. names or rownames have to be sample name. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase no genotype rsnp
#' @export
#' @return list to input to neg.beta.noGT.rsnp.priors.eff.stan
#' in.neg.beta.noGT.eff()

in.neg.beta.noGT.eff <- function(x,covar=1){
    
    N <- length(x$NB$counts) # number of individuals for NB    
    cov <- stan.cov(N,covar)  ## format covariates for stan
    ## g.NB
    g.NB <- lapply(x$NB$p.g, names)
    g.NB <- as.numeric(do.call(c,g.NB))

    ## p.NB
    p.NB <-  unname(do.call(c,unname(x$NB$p.g)))

    ## s.NB number of genotypes for each sample
    s.NB <- unname(sapply(x$NB$p.g,length))


    # n. counts    
    n.v <- do.call(c,unname(x$ase$n))
    s <- unname(sapply(x$ase$n,length))
    L <- sum(s)
    #p(H)
    p.v <- do.call(c,x$ase$p)
    ##gase
    gase <- unname(do.call(c,unname(x$ase$g)))
    # v
    v <- 0: max(x$ase$m)

    ## Index individuals according to having NB info and ASE
    ASEi <- as.numeric(names(Y) %in% names(x$ase$m))
    
    LL=list(N=N, G=sum(s.NB), A=length(x$ase$m), L=L, K=ncol(cov)-1, M=length(v), Y=Y ,sNB=s.NB, gNB=g.NB, pNB=p.NB, gase= gase, m=x$ase$m, n=n.v, pH=p.v, s=s, v=v, cov=cov, ASEi=ASEi)

    return(LL)
}


#' prepare input for stan noGT.test.stan, only uses individuals with ASE counts for NB, better to use in.neg.beta.noGT.eff
#'
#' This function allows you to prepare inputs for stan stan noGT.test.stan
#' @param x output list from stan.full.no.gt
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase no genotype rsnp
#' @export
#' @return list to input to neg.beta.noGT.rsnp.priors.eff.stan
#' in.neg.beta.noGT.eff()

in.noGT.test.eff <- function(x,covar=1){
    
    N <- length(x$NB$counts) # number of individuals
    if(length(covar)==1){
        cov <- matrix(1,nrow=N, ncol=2) # I need at least 2 cols so stan will recognise covar as a matrix, which helps with downstream calculations in stan.
     
    } else { #first extra column will be ignored
        cov <- cbind(rep(1, N), covar)
    }
    ## g.NB
    g.NB <- lapply(x$NB$p.g, names)
    g.NB <- as.numeric(do.call(c,g.NB))

    ## p.NB
    p.NB <-  unname(do.call(c,unname(x$NB$p.g)))

    ## s.NB number of genotypes for each sample
    s.NB <- unname(sapply(x$NB$p.g,length))
    
    # n. counts    
    n.v <- do.call(c,unname(x$ase$n))
    s <- unname(sapply(x$ase$n,length))
    L <- sum(s)
    #p(H)
    p.v <- do.call(c,x$ase$p)
    ##gase
    gase <- unname(do.call(c,unname(x$ase$g)))
    # v
    v <- 0: max(x$ase$m)
    LL=list(N=N, G=sum(s.NB), L=L, K=ncol(cov)-1, M=length(v), Y=unlist(x$NB$counts),sNB=s.NB, gNB=g.NB, pNB=p.NB, gase= gase, m=x$ase$m, n=n.v, pH=p.v, s=s, v=v, cov=cov)
    return(LL)
}



#' function for formatting input for trecase stan no GT, wrapping p.hap.pair and stan.neg.beta.prob.eff to avoid creating big objects when returning from p.hap.pair
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual. Starts from refernce panel input, calculates p(H|G) and then prepares stan input for full model
#' @param counts.g data table 1 row with total gene counts of gene per sample
#' @param rp.1r numeric of haplotypes from reference panel for 1 rsnp
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param c.ase output from tot.ase for the rsnp id
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.n minimun number of individuals with sufficient ase counts, trease default 5
#' @keywords stan input
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.trecase.noGT.eff()

stan.trecase.noGT.eff <- function(counts.g, rp.1r, rp.f, f.ase, c.ase , min.ase=5, min.ase.n=5){
  
    rp2 <- t(rbind(rp.f,rp.1r)) ## get ref panel fsnps and 1 rsnp in the same format I have functions from simulations
    
    ## calculate P(H|G) ref panel for the rsnp                   
    tmp <- p.hap.pair.s(rp2) 
    
    ## get haps for fsnps 
    h.samp <- hap_sam(x=f.ase)

    if(nrow(h.samp[[1]])==1) {
        h.samp <- lapply(h.samp,t)
    }
    
    h.samp = lapply(h.samp, function(x) {row.names(x)= rownames(c.ase); return(x)})
    
    tmp <- stan.full.no.gt(counts.g , M=tmp, l=h.samp, m=c.ase, ase=min.ase, n=min.ase.n)
    
    return(tmp)
    
    
}

#' function for formatting input for trecase stan no GT and missing values, prepares inputs only dependent on fsnps to feed stan.trecase.rna.noGT.eff
#'
#' This function allows you to format input for stan related to fsnps only
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param c.ase output from tot.ase for the rsnp id
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.n minimun number of individuals with suffiicent ase counts
#' @param ai data table with allelic imbalance estimate for reference panel bias correction, defaults to NULL
#' @keywords stan input
#' @export
#' @return list of 1) m= numeric with total ase counts. 2) n= list of matrices: each matrix has n.counts, rows samples, cols haplotypes compatible with the reference panel. 3) NB= list of data tables, each data table summarizing samples with the haplotypes compatible with RP, to make NB input, 4) f.comb, list relating names on n list and NB list to fsnps used for making haplotypes
#'
#' fsnp.prep2()

fsnp.prep2 <- function(rp.f, f.ase, c.ase , min.ase=5, min.ase.n=5, ai=NULL){ 
    
    ##  group individuals with no missing GT for the same fsnps ####
    stan.f <- fsnp.prep(rp.f, f.ase, c.ase , min.ase, min.ase.n)
    
    if(is.character(stan.f)){return(stan.f)}
    
    mis.sample <- stan.f$mis.sample
    u.f <- stan.f$u.f
    f.comb <- stan.f$f.comb
    
    ##  index samples with GT info for same fsnps
    index <- apply(u.f,2, function(i) apply(mis.sample, 2, function(j) sum(j==i)==nrow(f.ase)))
    rownames(index) <-  gsub("_GT","", rownames(index))

    ## check if I have sufficient individuals with sufficient counts before proceeding

    m.col <- grep("\\.m",colnames(c.ase), value=T)
    m.counts <-  rowSums(c.ase[,m.col, drop=FALSE])
    
    ## select ASE input when total ase counts are above threshold
    A <- which(m.counts>=min.ase)
    if(length(A)<min.ase.n){ ## at least n individuals with sufficient ase counts
           return("Not enough individuals with ASE counts")
    } 

##### based on f.comb, u.f and index prepare inputs with no missing values for stan.fsnp.noGT.eff

    fsnp.in <- lapply(1:length(f.comb), function(i) {
        samp <- rownames(index[index[,i],,drop=FALSE])
        cols.fase <- unlist(lapply(samp, grep, x=names(f.ase), value=TRUE))
        fsnps <- f.comb[[i]]
        cols.case <- unlist(lapply(fsnps, grep, x=colnames(c.ase), value=TRUE))
        ## for each entry I ignore checking for ase, it is only a subset of individuals and it was checked globally above
        tmp <- stan.fsnp.noGT.eff(rp.f=rp.f[fsnps,,drop=FALSE], f.ase=f.ase[id %in% fsnps, cols.fase, with=FALSE], c.ase=c.ase[samp, cols.case,drop=FALSE], NB="yes",min.ase= min.ase, min.ase.n=0,  ai=ai)
        return(tmp)
    })
    names(fsnp.in) <- names(f.comb)
    
    ## remove character elements (GT not compatible with ref panel)
    w <- sapply(fsnp.in, is.character)

    if(all(w)){
        return(fsnp.in[[1]])
    }
    
    if(any(w)){
        fsnp.in <- fsnp.in[!w]        
    }

    ## prepare outputs

    ## m counts
    m.out <- Reduce(c,lapply(fsnp.in, `[[`, "m"))
    ## remove NULL, same inds may have NB but not m
    m.out[sapply(m.out, is.null)] <- NULL

    ## n counts, list names link to fsnp used by f.comb
    n.out <- lapply(fsnp.in, `[[`,"n")
    ## remove NULL, same inds may have NB but not m
    n.out[sapply(n.out, is.null)] <- NULL

    if(!is.null(ai)) {
        ai.out <- lapply(fsnp.in, `[[`,"ai")
        vai.out <- lapply(fsnp.in, `[[`,"vai")
        ## remove NULL, same inds may have NB but not m
        ai.out[sapply(ai.out, is.null)] <- NULL
        vai.out[sapply(vai.out, is.null)] <- NULL
    }
    
    ## NB side
    nb <- lapply(fsnp.in, `[[`,"NB")

   

 
    if(!is.null(ai)){
        return(list(m=m.out, n=n.out, ai=ai.out, vai=vai.out, NB=nb, f.comb=f.comb))

    } else {
        
        return(list(m=m.out, n=n.out, NB=nb, f.comb=f.comb))

    }

    
}

#' helper function for formatting input for trecase stan no GT and missing values, prepares inputs only dependent on fsnps to feed stan.trecase.rna.noGT.eff
#'
#' This function allows you to format input for stan related to fsnps only
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param c.ase output from tot.ase for the rsnp id
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.n minimun number of individuals with suffiicent ase counts
#' @keywords stan input
#' @export
#' @return list of 1) mis.sample, logical matrix, for each sample indicates which fnsps are missing (TRUE) 2) u.f, matrix with unique cols from mis.samples, 3) f.comb list to relate u.f with name of snps, each element is a vector with name colname in u.f, and values names of fsnps that are not missing.
#'
#' fsnp.prep()

fsnp.prep <- function(rp.f, f.ase, c.ase , min.ase=5, min.ase.n=5){ 
    
#######  group individuals with no missing GT for the same fsnps ####
    ## get GT info
    GT <- f.ase[,grep("_GT$",names(f.ase)),with=F]

    ## identify missing GT by sample, MISSING==TRUE
    mis.sample <- apply(GT,2, function(i) i %in% "./.")
    if(!is.matrix(mis.sample)){
        cols <- names(mis.sample)
        mis.sample <- matrix(mis.sample, nrow=1)
        colnames(mis.sample) <- cols
    }
    
    rownames(mis.sample) <- f.ase$id
  

    ## remove samples with all missing data
    mis.sample <- mis.sample[,colSums(mis.sample)< nrow(GT), drop=F]

    if(ncol(mis.sample) < min.ase.n) return(paste("STOP: genotype information in less than" ,min.ase.n , "samples"))

    ## get unique combinations of fsnps
    u.f <- t(unique(t(mis.sample)))
    
    ## get the names of fsnps for each combination in u.f
    f.comb <- apply(u.f,2,function(i) names(i[which(!i)]))

    ## make f.comb list when is matrix (only one configuration)

    if(is.matrix(f.comb)){
        f.comb <- list(c(f.comb))
        names(f.comb) <- colnames(u.f)
    }

    return(list(mis.sample=mis.sample, u.f=u.f, f.comb=f.comb))
}




#' function for formatting input for trecase stan no GT and missing values, wrapping p.hap.pair and stan.neg.beta.prob.eff to avoid creating big objects when returning from p.hap.pair
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual. Starts from refernce panel input, calculates p(H|G) and then prepares stan input for full model
#' @param counts.g data table 1 row with total gene counts of gene per sample
#' @param rp.1r numeric of haplotypes from reference panel for 1 rsnp
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param c.ase output from tot.ase for the rsnp id
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.n minimun number of individuals with suffiicent ase counts
#' @param stan.f list with information for fsnps per sample, output from snp.prep
#' @keywords stan input
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.trecase.rna.noGT.eff()

stan.trecase.rna.noGT.eff <- function(counts.g, rp.1r, rp.f, f.ase, c.ase , min.ase=5, min.ase.n=5, stan.f){ 
    
#######  group individuals with no missing GT for the same fsnps ####
    mis.sample <- stan.f$mis.sample
    u.f <- stan.f$u.f
    f.comb <- stan.f$f.comb
    
    rp.comb <- lapply(1:length(f.comb), function(i) t(rbind(rp.f[rownames(rp.f) %in% f.comb[[i]],,drop=F] ,rp.1r))) ## get ref panel fsnps and 1 rsnp in the same format I have functions from simulations
    
    ## calculate P(H|G) ref panel for the rsnp, according to available fsnps             
    tmp <- lapply(rp.comb, p.hap.pair.s)
    
    ## get haps for fsnps and rsnp, need to index samples with GT for same fsnps
    index <- apply(u.f,2, function(i) apply(mis.sample, 2, function(j) sum(j==i)==nrow(f.ase)))
    
    h.samp <- lapply(1:ncol(index), function(i) {
        x=f.ase[id %in% rownames(u.f)[(!u.f[,colnames(index)[i]])]  ,names(index[,i])[index[,i]],with=F]
        l=hap_sam(x)
        if(nrow(x)==nrow(l[[1]])){ ## l needs to be transposed relative to x
            l <- lapply(l,t)
        }
        return(l)
    })
    
    ##get sample names for samples used for each combination of fsnps to aid in preparing inputs
    samps <- lapply(1:ncol(index), function(i) gsub("_GT","", names(index[,i])[index[,i]]))
    
    h.samp = lapply(1:length(h.samp), function(x) lapply(h.samp[[x]], function(i) {row.names(i)= samps[[x]]; return(i)}))

    
    ## prepare inputs, use n=0 because was tested globally above, select samples and snps for each haplotype
    tmp1 <- lapply(1:length(h.samp), function(i) stan.full.no.gt(counts=counts.g[,samps[[i]],with=F], M=tmp[[i]], l=h.samp[[i]], m=c.ase[samps[[i]], unlist(lapply(f.comb[[i]], function(j) grep(j, colnames(c.ase)))) ,drop=F], ase=min.ase, n=0,rna="yes"))
    
    if(length(tmp1)==1){ ## remove error message return for all samples
        return(tmp1)
    }
    ## Remove elements from tmp1 with mo information
    tmp1 <- tmp1[sapply(tmp1,function(i) class(i)!="character")]
    
    ## combine inputs for stan into one list, group NB side and ASE side
    ## NB side
    counts.all <- Reduce(cbind,lapply(tmp1, function(x) x$NB$counts))
    p.g.all <- Reduce(c, lapply(tmp1, function(x) x$NB$p.g))
    
    ## ASE side, remove elements from tmp1 with no ase info
    tm1 <- tmp1[sapply(tmp1,function(i) class(i$ase)!="character")]
    ase.all <- lapply(1:4, function(i) Reduce(c, lapply(tm1, function(x) x$ase[[i]])))
    names(ase.all) <- c("m","n","p","g")                  
    tmp2 <- list(NB=list(counts=counts.all, p.g= p.g.all), ase=ase.all)
    return(tmp2)
    
    
}

#' function for formatting input for trecase stan no GT and missing values, wrapping p.hap.pair and stan.neg.beta.prob.eff to avoid creating big objects when returning from p.hap.pair, v2
#'
#' This function allows you to format input for stan per individual: constructs list of vectors with counts/ase/p(H|G) for all haplotypes concordant with "g", each vector one individual. Starts from refernce panel input, calculates p(H|G) and then prepares stan input for full model
#' @param counts.g data table 1 row with total gene counts of gene per sample
#' @param rp.1r numeric of haplotypes from reference panel for 1 rsnp
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param c.ase output from tot.ase for the rsnp id
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.n minimun number of individuals with suffiicent ase counts
#' @param stan.f list with information for fsnps per sample, output from snp.prep
#' @keywords stan input
#' @export
#' @return list of 1) data table with total counts (y) and genotype rsnp (g) 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n. The first haplotype pair corresponds to {h1,h2}, the "true" one.
#'
#' stan.trecase.rna.noGT.eff2()

stan.trecase.rna.noGT.eff2 <- function(counts.g, rp.1r, rp.f, stan.f) { 
    
#######  group individuals with no missing GT for the same fsnps ####
    m <- stan.f$m
    n <- stan.f$n
    NB <- stan.f$NB
    f.comb <- stan.f$f.comb
    if("ai" %in% names(stan.f)) {
        ai0 <- stan.f$ai
        vai0 <- stan.f$vai
    }
    

    ## Prepare inputs for stan.rsnp.noGT.eff by fsnps
    nb <- list()
    n.list <- list()
    p.list <- list()
    g.list <- list()
    ai0.list <- list()
    vai0.list <- list()
    
    for( i in 1:length(f.comb)) {
        NB.sub <- NB[names(NB)==names(f.comb)[i]]
        n.sub <- n[names(n)==names(f.comb)[i]]
        if("ai" %in% names(stan.f)){
            ai0.sub <- ai0[names(ai0)==names(f.comb)[i]]
            vai0.sub <- vai0[names(vai0)==names(f.comb)[i]]
        }
        
        if(length(NB.sub) | length(n.sub)){        
            rp.f2 <- rp.f[f.comb[[i]],]
            rp2 <- t(rbind(rp.f2,rp.1r))
            M <- p.hap.pair.s(rp2)
            M.cond <- mat.col(M)
        }
        
        if(length(NB.sub)){ 
            nb[[i]] <- stan.nb(counts.g,M,M.cond,NB=NB.sub[[1]])
          
        }
        if(length(n.sub)){
            ## make named vector with names rownames(M) and value swapped rownames(M) once
            row.split <- unlist(strsplit(rownames(M), ","))
            row.2 <- row.split[1:length(row.split) %%2 == 0]  ## take even entries
            row.1 <- row.split[1:length(row.split) %%2 ==1] ## odd entries
            s.M <- paste(row.2, row.1, sep=",")  ## invert haplotypes
            names(s.M) <-  rownames(M)
            
            if("ai" %in% names(stan.f)){
                if(length(ai0.sub)){
                    tmp <- Reduce(c, mapply(function(a,b,c) stan.ase(M, M.cond,s.M, n.mat=a, ai0.mat=b, vai0.mat=c),
                                            a=n.sub[[1]],
                                            b=ai0.sub[[1]],
                                            c=vai0.sub[[1]]))  ## simplify

                    n.list <- c(n.list, tmp[grep("n", names(tmp))])
                    p.list <- c(p.list, tmp[grep("p", names(tmp))])
                    g.list <- c(g.list, tmp[grep("g", names(tmp))])                            
                    ai0.list <- c(ai0.list, tmp[grep("ai0", names(tmp))])
                    vai0.list <- c(vai0.list, tmp[grep("vari0", names(tmp))])
                
                }
            } else {
                tmp <- Reduce(c, lapply(n.sub[[1]], function(j) stan.ase(M, M.cond,s.M, n.mat=j)))  ## simplify
                n.list <- c(n.list, Reduce(c, lapply(tmp, `[`, "n")))
                p.list <- c(p.list, Reduce(c, lapply(tmp, `[`, "p")))
                g.list <- c(g.list, Reduce(c, lapply(tmp, `[`, "g")))
                
            }
            
            
        }
             
        
    }
    
    ## Arrange lists to return ##
    
    ## NB remove null elements
    nb[sapply(nb, is.null)] <-NULL
    tot.counts <- Reduce(cbind,lapply(nb, function(i) i$counts))
    p.g <- Reduce(c, lapply(nb, function(i) i$p.g))
    names(p.g) <- names(tot.counts)

    NB <- list(counts=tot.counts, p.g=p.g)

    ## ASE ##
    n.list <- Reduce(c,n.list)
    p.list <- Reduce(c, p.list) 
    g.list <- Reduce(c, g.list)

    
    ASE <- list(m=m, n=n.list[names(m)], p=p.list[names(m)], g=g.list[names(m)])

    if("ai" %in% names(stan.f)) {
        ai0.list <- Reduce(c, ai0.list)
        ASE$ai0 <- ai0.list[names(m)]
        vai0.list <- Reduce(c, vai0.list)
        ASE$vai0 <- vai0.list[names(m)]
    }
    

    return(list(NB=NB, ase=ASE))

}

    

        
    
                   
#' Fix GT for rsnp on negative binomial input only for QC purposes in stan neg.beta.noGT.rsnp.priors.eff.stan
#'
#' This function allows you to fix the genotype of the rsnp for stan neg.beta.noGT.rsnp.priors.eff.stan input, only for the negtive binomial side for QC purposes
#' @param x input for neg.beta.noGT.rsnp.priors.eff.stan, one rsnp
#' @param y data table with genotypes for one rsnp, may be coded as 0,1,-1,2, output from rec_mytrecase_rSNPs. Cols inlcude those named as sample_GT, plust others describing the snp 
#' @keywords stan fix genotype no genotype rsnp
#' @export
#' @return list to input to neg.beta.noGT.rsnp.priors.eff.stan with element $pg modified to fix genotype p=1
#' fix.noGT()

fix.noGT <- function(x,y){
    ## select samples in x present in y
    y.sub <- y[,which(names(y) %in% paste0(names(x$Y), "_GT")), with=F]
    N=length(x$Y)
    ## gNB is now y.sub, all the rest is the same
    LL=list(N=N, G=N, L=x$L, K=x$K, M=x$M, Y=x$Y, sNB=rep(1, N), gNB=unname(unlist(y.sub)), pNB=rep(1,N), gase=x$gase, m=x$m, n=x$n, pH=x$pH, s=x$s, v=x$v, cov=x$cov)
    return(LL)
}

                       
#' prepare input for stan neg.beta.noGT.rsnp.priors.eff2.stan
#'
#' This function allows you to prepare inputs for stan neg.beta.noGT.rsnp.priors.eff2.stan
#' @param x output list from stan.full.no.gt
#' @param covar matrix or numeric with covariates, rows individuals, cols covariates. Individuals in same order as in x, do not include col of ones for intercept. Defaults to intercept only
#' @keywords stan input trecase no genotype rsnp
#' @export
#' @return list to input to neg.beta.noGT.rsnp.priors.eff.stan
#' in.neg.beta.noGT.eff2()

in.neg.beta.noGT.eff2 <- function(x,covar=1){
    
    N <- length(x$NB$counts) # number of individuals for NB

    ## Count haplotypes per genotype of cis-SNP

    gNB.l <- lapply(x$NB$p.g, function(i) as.numeric(names(i)))
    

    ## select elements in gNB.l with ASE info
   gNB.ase  <- gNB.l[names(x$ase$g)]

    ## sort each indvidual ASE info (gase, n , p, ai0) in order of gNB, so to have the haps ordered by genotype

    ## get indexes for sorting based on g 
    indx <- mapply(function(a,b) order(match(abs(b), a)),
                    a=gNB.ase, b=x$ase$g, SIMPLIFY=F)                  
     
    ## apply indx to relevant ase cols
    
    gase <- aux.sort(x$ase$g, indx)
    n.v  <- aux.sort(x$ase$n, indx)
    p.v <- aux.sort(x$ase$p, indx) 
             
    if(any(names(x$ase) == "ai0")) {
        ai0 <- aux.sort(x$ase$ai0, indx)
        vai0 <- unlist(unname(x$ase$vai0))  ## doesnt need sorting, unique value per individual
    }
    
   
    ##for those inds with ASE get the number of haps for each genotype, I dont need them in order, I am just counting so I can use x$ase$inputs
    h2g <- mapply(function(a,b) {
        table(abs(b)[abs(b) %in% a])
    }
 ,  a=gNB.ase, b=x$ase$g, SIMPLIFY=F)

    ## for those without ASE add 0 to each genotype, need to make "0" as otherwise dropped to one value when unlisting below

    h2g.no <- lapply(gNB.l[!names(gNB.l) %in% names(x$ase$g)], function(i) {
        v=rep("0", length(i))
        names(v)=i
        return(v)
        })

    h2g <- c(h2g, h2g.no)

    ## need to sort NB with first ase names  and then no ase names(as in names(h2g)) ##

    ## includes, g, p , counts and covar
    g.NB <- as.numeric(unlist(lapply(x$NB$p.g[names(h2g)] , names))) ## genotypes are in names, first order elements by ase info
    p.NB <-  unname(unlist(x$NB$p.g[names(h2g)]))
    Y=unlist(x$NB$counts[,names(h2g), with=F])
    if(is.matrix(covar)){  ## make sure covar has the same samples as x$NB$counts
        covar = covar[names(Y),]
    }

    ## Index individuals according to having NB info and ASE, first col 1 if NB and ASE, second col indexes ASE individuals
    ASEi <- as.numeric(names(Y) %in% names(x$ase$m))
    mat <- matrix(c(ASEi, cumsum(ASEi)), ncol=2)
    mat[mat[,1]==0,2] <- 0

    
     ## prepare stan inputs
    s <- unname(sapply(x$ase$n,length))
    L <- sum(s)   
    cov <- stan.cov(N,covar)  ## format covariates for stan
    ## s.NB number of genotypes for each NB sample
    s.NB <- sapply(x$NB$p.g[names(h2g)],length)


    ## make it unnamed vector to return
    h2g <- as.numeric(unlist(h2g))

   
    LL=list(N=N, G=sum(s.NB), A=length(x$ase$m), L=L, K=ncol(cov)-1, Y=Y,sNB=s.NB, gNB=g.NB, pNB=p.NB, gase= gase, m=x$ase$m, n=n.v, pH=p.v, s=s, cov=cov, ASEi=mat, h2g=h2g)


    if(any(names(x$ase) == "ai0")){
        ## convert vai0 to sdai0 as stan takes standard deviation as input for normal distribution
        LL <- c(LL, list(ai0=ai0, sdai0=sqrt(vai0)))

    }
        
    return(LL)
}

#' Aux function to sort ase according by ascending genotype order (based on NB genotype)
#'
#' @param a list to sort
#' @param ind list with indexes
#' @keywords sort ase by Gi=g
#' @export
#' @return vector to input for stan
#'
#' aux.sort()

aux.sort <- function(a,b){
    tmp <- unlist(mapply(function(x,y) x[y], x=a, y=b, SIMPLIFY=F, USE.NAMES=F))
    return(tmp)
    }
    



###########################################################################################################################################################################
######## Function to avoid running issues with many models in stan

#' deal with too many ddls loaded when running stan many times: https://github.com/stan-dev/rstan/issues/448
#'
#' This function allows you to unload not in use ddls
#' @param model name of object returned by stan_model, to avoid recompilation.
#' @keywords stan bug ddls
#' @export
#' @return unload old ddls
#' unload.ddl()

unload.ddl <- function(model){
    dso_filename = model@dso@dso_filename
    loaded_dlls = getLoadedDLLs()
    if (dso_filename %in% names(loaded_dlls)) {
        ##message("Unloading DLL for model dso ", dso_filename)
        model.dll = loaded_dlls[[dso_filename]][['path']]
        dyn.unload(model.dll)
    } else {
        ##message("No loaded DLL for model dso ", dso_filename)
    }
                    
    loaded_dlls = getLoadedDLLs()
    loaded_dlls <- loaded_dlls[grep("^file", names(loaded_dlls),value=T)]
    if (length(loaded_dlls) > 10) {
        for (dll in head(loaded_dlls, -10)) {
            ##message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
            dyn.unload(dll[['path']])
        }
    }
}


################################################################################################################
#################### New input formulation for stan #############################################

#' function for formatting input for trecase stan no GT, deals with haps, n counts and m counts of fSNPs to streamline computation, optional to include allelic imbalance estimate for dealing with reference panel bias
#'
#' This function allows you to format input for stan per individual: precomputes m counts and n counts for each individual based on haplotypes compatible with genotypes of the fsnps based on the reference panel
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param f.ase data.table with GT and ASE for fSNPs, output from gt.as, same fsnps in same order as rp.f
#' @param c.ase output from tot.ase for the rsnp id
#' @param min.ase ase cut-off default is 5 counts, trecase default 5
#' @param min.ase.n minimun number of individuals with sufficient ase counts, trease default 5
#' @param NB, whether to return NB element of output, defalts to "yes" as it is required for noGT, use "no" for GT input preparation
#' @param ai data table with AI estimates per SNP, defaults to NULL
#' @keywords stan input fsnp
#' @export
#' @return list of 1) m= numeric with total ase counts. 2) n= list of matrices: each matrix has n.counts, rows samples, cols haplotypes compatible with the reference panel. 3) NB= data table summarizing samples with the haplotypes compatible with RP, to make NB input, 
#'
#' stan.fsnp.noGT.eff()

stan.fsnp.noGT.eff <- function(rp.f, f.ase, c.ase , NB="yes", min.ase=5, min.ase.n=5, ai=NULL){

    ## get haps for fsnps 
    h.samp <- hap_sam2(x=f.ase, w=rownames(c.ase))
    
    ## get GT for fsnps
    g <- Reduce("+",h.samp)
    g <-apply(g[,,drop=F],1,paste,collapse="")  
    ## Select samples with GT compatible with reference panel

    ## calculate P(H|G) ref panel for the fsnps,
    M <- p.hap.pair.s(t(rp.f))

    g.all.comp <- g[g %in% colnames(M)]
    if(length(g.all.comp) == 0 ){
        return("GT fsnps is not compatible with Ref panel")
    }
    ## NB side: ALL INDIVIDUALS with GT compatible with RP, get haps
    if(NB=="yes"){
        if(nchar(g.all.comp[1])==1){ ## one fsnp, tested in first element of g.all.comp
            h.short.all <- lapply(h.samp, function(i) i[names(g.all.comp),])
        } else {
            h.short.all <- lapply(h.samp, function(i) apply(i[names(g.all.comp),,drop=FALSE],1,paste,collapse=""))
        }
        
        ## summary by unique(haplotype) indicating in which samples are found and in the row of M that the hap pair is found
        u.haps.all <- sum.g(h.short.all)
        ## recode to sample names instead of M row

        
    }
    

    ## ASE side: samples with sufficient ASE at least in n individuals and compatible with ref panel
    
    ## get GT for fsnps in samples with sufficient ASE at least in n individuals and compatible with ref panel
    m.col <- grep("\\.m",colnames(c.ase), value=T)

    
    m.counts <-  rowSums(c.ase[,m.col, drop=FALSE])
    

    ## keep only samples with GT compatible with RP

    m.counts <- m.counts[names(g.all.comp)]
    
    ## select ASE input when total ase counts are above threshold
    A <- which(m.counts>=min.ase)
    
    if(length(A)<min.ase.n){ ## at least n individuals with sufficient ase counts (unknown genotype)
        if(NB=="yes"){
            return(list(NB=u.haps.all))
        } else {
            return("Not enough individuals with ASE counts")
        }
        
    }

    g.comp <- g[names(A)]
    #g.obs <-apply(g[A,,drop=F],1,paste,collapse="")    
    ##names(g.obs)=names(A)
    ## select inds that have GT fsnps compatible with ref panel
    ##g.comp <- g.obs[g.obs %in% colnames(M)]

    #print(g.comp)
    if(!length(g.comp) ){
         if(NB=="yes"){
            return(list(NB=u.haps.all))
        } else {
            return("GT fsnps is not compatible with Ref panel")
        }
         
    }
    ## proceed to calculate m.counts: only inds with ASE info and GT compatible with ref panel
    m.trim <- c.ase[names(g.comp),,drop=FALSE] 
    if(nchar(g.comp[1])==1){ ## one fsnp, tested in first element of g.comp
        h.short <- lapply(h.samp, function(i) i[names(g.comp),])
        m.counts <- m.trim[,m.col]
        if(!length(names(m.counts))){ names(m.counts) <- names(g.comp)}
    } else {
        h.short <- lapply(h.samp, function(i) apply(i[names(g.comp),,drop=FALSE],1,paste,collapse=""))
        m.counts <-  rowSums(m.trim[,m.col,drop=FALSE])
    }

    ## Prepare list with m.counts, n.counts and the corresponding hap pair
    inp <- list(m=m.counts,n=list())

    ## Prepare list with ai estimate and vector for variance ai estimate
    if(!is.null(ai)) {
        inp$ai = list()
        inp$vai = list()
        ## add logit column
        ai[, logitAI_post:=qlogis(AI_post)]
    }
    

    ## to make it faster: get inputs for unique(haplotypes) and then expand it to all samples
    u.haps.dt <-sum.haps(h.short)

    ## make hap1 and hap2 character columns if they arent

    u.haps.dt[ , hap1:=as.character(hap1)][,hap2:= as.character(hap2)]
        
    ## get indexes for u.haps in hap.dt
    if(nrow(u.haps.dt)==1){## general form doesnt work only accepts one value
        idx=as.numeric(unlist(strsplit(u.haps.dt$ind, ",")))
    } 
    ## add hap1.2 col to u.haps.dt
    u.haps.dt[, hap1.2:=paste(hap1,hap2,sep=",")]
    ## add row number that corresponds to hap1.2 in M, if misssing row=0
    u.haps.dt[, row.M:=hap.p.no.gt(M,hap1.2)]

    
    for(i in 1:nrow(u.haps.dt)) {
        if(nrow(u.haps.dt)>1){
            idx <- unlist(u.haps.dt$ind[i]) ## index individuals with same hap combinations
        }
        ## get matrix for n.counts, for all individuals with same haps, each row is the n for each individual, cols will be all possible haps.
        n.v.mat <- matrix(,nrow=length(idx), ncol=0)  ## for all individuals with same haps, each row is the n for each individual, cols will be all possible haps.
        n.col <- grep("\\.n",colnames(c.ase), value=T)
        
        ## get vector for  AI.null, for all individuals with same haps, values correspond to all possible haps.
        if(!is.null(ai)){
            ai.v.mat <- matrix(,nrow=length(idx), ncol=0)
            vai.mat <- matrix(,nrow=length(idx), ncol=0)
            ## get positions that AI needs to be swapped
            ## get positions with g=1
            hets <- unlist(gregexpr("1", u.haps.dt[i,g.fsnps])) ## first element of the list are the indices
        }
        
        col <- gsub("\\.m", "", m.col)
        if(u.haps.dt$row.M[i] >0){ ## obs hap in ref panel   
            n.v.mat <- cbind(n.v.mat, rowSums(m.trim[idx,n.col,drop=FALSE]))
            colnames(n.v.mat) <- u.haps.dt$hap1.2[i]
            ex <- u.haps.dt$row.M[i] ## haplotype to exclude from further analysis
            if(!is.null(ai)){ 
                ## check whether alt allele for the het is in hap1, if so AI is -AI (logist scale)
                swap <- unlist(strsplit(u.haps.dt[i,hap1], ""))[hets] == 1
                temp.ai <- ai[match(col,id), ]
                ## swap logist when swapping hap
                temp.ai[hets[swap], logitAI_post := -logitAI_post]
                temp2.ai <- temp.ai[,.(logitAI_post)]
                ai.v.mat <- cbind(ai.v.mat, m.trim[idx,m.col,drop=FALSE] %*% as.matrix(temp2.ai)/rowSums(m.trim[idx,m.col,drop=FALSE]))
                colnames(ai.v.mat) <- u.haps.dt$hap1.2[i]

                ## add variance
                vai.mat <- cbind(vai.mat,   m.trim[idx,m.col,drop=FALSE]^2 %*% (1/temp.ai$NREF_post + 1/temp.ai$NALT_post)/rowSums(m.trim[idx,m.col,drop=FALSE])^2 )
                colnames(vai.mat) <- u.haps.dt$hap1.2[i]
                
            }
            
        }       
        ## get all other possible haplotypes compatible with genotype of fsnps from ref panel
        if(exists("ex")){
            hap.others <- names(which( M[-ex,which(colnames(M) %in% u.haps.dt$g.fsnps[i])]>0))
        } else {
            hap.others <- names(which( M[,which(colnames(M) %in% u.haps.dt$g.fsnps[i])]>0))
        }
        ## check if they are any
        if(length(hap.others)){
            hap1 <- u.haps.dt$hap1[i]
            hap.x <- lapply(hap.others, function(i) unlist(strsplit(i,","))[1])  ## hap 1
            for(k in 1:length(hap.x)){
                ##get the differences between hapx hap1 to swap ase
                dif <- which(unlist(strsplit(hap1,split=""))!=unlist(strsplit(hap.x[[k]],split="")))
                ## get ase
                n.v.mat <-cbind(n.v.mat, rowSums(m.trim[idx,c(n.col[-dif],m.col[dif]), drop=FALSE]) - rowSums(m.trim[idx,n.col[dif], drop=FALSE]))

                if(!is.null(ai)){
                    ## check whether alt allele for the het is in hap1, if so AI is -AI (logit scale)
                    ## for this part hets remains the same, as the genotype is the same
                    swap <- unlist(strsplit(hap.x[[k]], ""))[hets] == 1
                    ## reset                   
                    temp.ai <- ai[match(col,id), ]
                    ## replace
                    temp.ai[hets[swap], logitAI_post := -logitAI_post]
                    temp2.ai <- temp.ai[,.(logitAI_post)]
                    ## temp.v <- sum(temp.ai$Total_post * temp.ai$lAI_post)/sum(temp.ai$Total_post)             
                    temp.mat <- m.trim[idx,m.col,drop=FALSE] %*% as.matrix(temp2.ai)/rowSums(m.trim[idx,m.col,drop=FALSE])
                    ## temp.mat <- matrix(temp.v, nrow=length(idx),dimnames=list(rownames( m.trim[idx,])))
                    colnames(temp.mat) <- NULL
                    ai.v.mat <- cbind(ai.v.mat, temp.mat)
                    vai.mat <- cbind(vai.mat,   m.trim[idx,m.col,drop=FALSE]^2 %*% (1/temp.ai$NREF_post + 1/temp.ai$NALT_post)/rowSums(m.trim[idx,m.col,drop=FALSE])^2 )

                    
                }
                

            }
        }
        colnames(n.v.mat) <- c(colnames(n.v.mat)[1:(ncol(n.v.mat)-length(hap.others))],hap.others)
        inp$n[[i]] <- n.v.mat

        if(!is.null(ai)){
            colnames(ai.v.mat) <- c(colnames(ai.v.mat)[1:(ncol(ai.v.mat)-length(hap.others))],hap.others)
            inp$ai[[i]]  <- ai.v.mat
            colnames(vai.mat) <- colnames(ai.v.mat)
            inp$vai[[i]] <- vai.mat
        }
        
    }

    
    if(NB=="yes"){
            inp$NB <- u.haps.all
    }
    ## some matrices in inp$n may have same col names but in different order, I will merge them as it wont affect assigning haplotype probabilities when adding rSNP

    ## only applies if more than 2 cols
    len <- lapply(inp$n, ncol)
    max.len <- max(unlist(len))

    if(max.len >1){       
        n.names <- lapply(inp$n, colnames)

        ## compare elements same number of cols, from 2 onwards
        ## when I can merge them I merge them
        ## keep record of merged matrices so I can delete them at the end so indices are not affected during the loop
        index <- c()
        ##colnames(index) <- c("j", "match")
        for(i in 2:max.len){
            w <- which(len==i)
            ## pair wise comparisons
            ## get j when cols have the same elements regardless order
            
            for (j in w){
                l <- lapply(sapply(n.names[w[-c(1:which(w==j))]], FUN=intersect, n.names[[j]]) , length) == i
                if(any(l)){                
                    ## get relevant matrices
                    inx <- c(j ,w[-c(1:which(w==j))][l])
                    ## sort cols same as the first one
                    for(k in 2:length(inx)){
                        inp$n[[inx[k]]] <- inp$n[[inx[k]]][,n.names[[j]], drop=F]
                        ## rbind
                        inp$n[[j]] <- rbind(inp$n[[j]], inp$n[[inx[k]]])
                        
                        if(!is.null(ai)) {  # same for ai
                            inp$ai[[inx[k]]] <- inp$ai[[inx[k]]][,n.names[[j]], drop=F]
                            inp$vai[[inx[k]]] <- inp$vai[[inx[k]]][,n.names[[j]], drop=F]
                            ## rbind
                            inp$ai[[j]] <- rbind(inp$ai[[j]], inp$ai[[inx[k]]])
                            inp$vai[[j]] <- rbind(inp$vai[[j]], inp$vai[[inx[k]]])
                        }
                    }
                    
                    ## sort rownames
                    inp$n[[j]] <- inp$n[[j]][names(m.counts)[which(names(m.counts) %in% rownames(inp$n[[j]]))], ]
                    index <- c(index, inx[2:length(inx)])

                    if(!is.null(ai)){
                        inp$ai[[j]] <- inp$ai[[j]][names(m.counts)[which(names(m.counts) %in% rownames(inp$ai[[j]]))], ]
                        inp$vai[[j]] <- inp$vai[[j]][names(m.counts)[which(names(m.counts) %in% rownames(inp$vai[[j]]))], ]
                    }
                    
                }
            }
        }
            
        ## remove redundant matrices
        inp$n[index] <- NULL
        if(!is.null(ai)) {
            inp$ai[index] <- NULL
            inp$vai[index] <- NULL
        }
        
        
    }
    
    return(inp)
}



#' helper function to summarise the haplotypes of fnsps across individuals
#'
#' This function allows you to format summarise the haplotypes of fnsps across individuals
#' @param h.short list with 2 elements, first character named vector with haplotype 1 for fsnps, if 3 fsps each element will look as "010", etc. names are sample ID. Second element corresponds to hap2
#' @keywords summary fsnps
#' @export
#' @return data table with columns: hap1, hap2, g.fsnps (GT fsnp), ind: row number from h.short for the individuals that have that hap combination, hap1.2, hap1 and hap2 separated by ",", row.M row number in which that hap combination is found in input M
#'
#' sum.haps()

sum.haps <- function(h.short){
 ## summary by unique(haplotype) indicating in which samples are found
    hap.s <- Reduce(cbind,h.short)
    colnames(hap.s) <- c("hap1", "hap2")
    hap.dt <- data.table(hap.s,keep.rownames=T)
    u.haps <- unique(hap.s,MARGIN=1)
    u.haps.dt <- data.table(u.haps)
    ## add geno fsnps
    if(nchar(h.short[[1]][1])==1) { ## 1 fsnp
        u.haps.dt[,g.fsnps:=hap1+hap2]
    } else {
        u.haps.dt[,g.fsnps:=Map(add.geno,hap1,hap2)]
    }
    
    ## get indexes for u.haps in hap.dt
    if(nrow(u.haps.dt)==1){## general form doesnt work only accepts one value
        u.haps.dt[, ind:= paste(which(hap.dt$hap1==hap1 & hap.dt$hap2==hap2), collapse=",")]        
    } else {            
        u.haps.dt[, ind:=lapply(1:nrow(u.haps.dt), function(i) which(hap.dt$hap1==hap1[i] & hap.dt$hap2==hap2[i]))]
    }
    return(u.haps.dt)
    
    ## add hap1.2 col to u.haps.dt
    #u.haps.dt[, hap1.2:=paste(hap1,hap2,sep=",")]
    ## add row number that corresponds to hap1.2 in M, if misssing row=0
    #u.haps.dt[, row.M:=hap.p.no.gt(M,hap1.2)]
   
}

#' helper function to summarise genotypes of fnsps across individuals
#'
#' This function allows you to format summarise the genotypes of fnsps across individuals
#' @param h.short list with 2 elements, first character named vector with haplotype 1 for fsnps, if 3 fsps each element will look as "010", etc. names are sample ID. Second element corresponds to hap2
#' @keywords summary fsnps
#' @export
#' @return data table with columns: hap1, hap2, g.fsnps (GT fsnp), ind: row number from h.short for the individuals that have that hap combination, hap1.2, hap1 and hap2 separated by ",", row.M row number in which that hap combination is found in input M
#'
#' sum.g()

sum.g <- function(h.short){
    ## summary by unique(haplotype) indicating in which samples are found
    if(is.numeric(h.short[[1]])){
        g.s <- Reduce("+",h.short)
    } else {
         g.s <- Reduce(add.geno,h.short)
    }
    g.m <- matrix(g.s, ncol=1, dimnames=list(names(g.s), "g.fsnps"))
    g.dt <- data.table(g.m, keep.rownames=T)
    u.g <- unique(g.m,MARGIN=1)
    u.g.dt <- data.table(u.g)
    
    ## get indexes for u.haps in hap.dt
    if(nrow(u.g.dt)==1){## general form doesnt work only accepts one value
        u.g.dt[, ind:= paste(rownames(g.m)[which(g.dt$g.fsnps==g.fsnps)], collapse=",")]        
    } else {            
        u.g.dt[, ind:=lapply(1:nrow(u.g.dt), function(i)rownames(g.m)[ which(g.dt$g.fsnps==g.fsnps[i])])]
    }
    return(u.g.dt)
    
   
}


#' function for formatting input for trecase stan no GT, deals with genotype of rsnp and p(Grsnp|G fsnp) & p(H|G) to make stan input
#'
#' This function allows you to format input for stan per individual: using m and n counts for each individual it adds the genotype of rsnp to prepare stan input
#' @param counts.g DT of total gene counts, colnames sample names
#' @param rp.f matrix of haplotypes from reference panel for fsnps, rows=snps, cols samples from reference panel
#' @param rp.1r matrix of phased GT for rsnp from reference panel, same format as for rp.f
#' @param stan.f list output of stan.fsnp.noGT.eff with m counts, n counts and p(H|G)
#' @keywords stan input rsnp
#' @export
#' @return list of stan input: NB: list of 1) data table with total counts (y) and list p.g, each element a named vector with p(genotypes) rsnp, names 0,1,2. ASE: 2) data.table with data for ASE: genotype rsnp (g) and total ase counts (m). 3) list of vectors n: each vector with ase n.counts compatible with each haplotype pair for 1 individual. 4) list of vectors p: each vector with p for each haplotype corresponding to n.
#'
#' stan.rsnp.noGT.eff()

stan.rsnp.noGT.eff <- function(counts.g, rp.f, rp.1r, stan.f) {
 
    rp2 <- t(rbind(rp.f,rp.1r)) ## get ref panel fsnps and 1 rsnp in the same format I have functions from simulations
    
    ## calculate P(H)ref panel for the fsnps+rsnp                   
    M.all <- p.hap.pair.s(h=rp2) 
    ## get input for NB, consider all individuals irrespective of ASE counts
    ## get matrix p(H|GT fsnps)
    M.all.cond <- mat.col(M.all)

    ## make data table with no 0 entries in M.all to get genotype of rsnp for a given hap.pair
    gt.rf <- data.table(summary(M.all))
    setkey(gt.rf,i)    

    ## for each stan.f$NB$g.fsnps get p(gf+rsnp| gfsnp)
    u <- unique(as.character(stan.f$NB$g.fsnps))
    p.g <- lapply(u,  function(i) {
        w <-  which(M.all.cond[,i]!=0)
        tmp <- M.all.cond[w,i]
        ## for each hap pair in names(tmp) get their corresponding genotype (colname with no-0 entry)
        gt.rf2 <- colnames(M.all)[gt.rf[i %in% w,j]]
        
        names(tmp) <- substr(gt.rf2,nchar(gt.rf2[1]), nchar(gt.rf2[1])) ## last character from gt.rf is geno rsnp
        tmp <- tapply(tmp, names(tmp), sum) ## add elements in vector with same genotype
        return(tmp)
    })
    names(p.g) <-u

    
    ## assign p.g to each sample
    u.haps.all <- stan.f$NB
    if(nrow(stan.f$NB)==1){## only one genotype
        
        samp.ind <- sort(unlist(strsplit(stan.f$NB$ind, ",")))
        
        p.g.comp <- rep(list(p.g[[1]]), length(samp.ind))
    } else {
        samp.ind <- sort(unlist(u.haps.all[, ind]))
        p.g.comp <- rep(list(p.g[[1]]), length(samp.ind)) ## start with first element and replace
        for(i in seq_along(names(p.g))) {
            ind <- sort(unlist(stan.f$NB[g.fsnps==names(p.g)[i],ind]))
            p.g.comp[ind] <- rep(list(p.g[[i]]),length(ind))
        }
    }
    
    ##  total counts for relevant samples
    ## names in counts include gene_id, need to shift samp.ind by 1
    y <- counts.g[,(samp.ind+1),with=F]

    
    ##ASE
    n.all <- list()
    g.all <- list()
    pH.all <- list()
    ## make named vector with names rownames(M) and value swapped rownames(M)
    row.split <- unlist(strsplit(rownames(M.all), ","))
    row.2 <- row.split[1:length(row.split) %%2 == 0]  ## take even entries
    row.1 <- row.split[1:length(row.split) %%2 ==1] ## odd entries
    s.M <- paste(row.2, row.1, sep=",")  ## invert haplotypes
    names(s.M) <-  rownames(M.all)

    for(i in 1:length(stan.f$n)) {
        n.mat <- stan.f$n[[i]]
        ## get hap pairs for fsnps and potentially any rsnp
        hap.f.p <- lapply(colnames(n.mat), function(i) unlist(strsplit(i,",")))
        hap.r.p <- lapply(hap.f.p, function(i) sapply(i, function(j) paste0(j,0:1)))
        hap.r.p <- lapply(hap.r.p, function(i) as.vector(outer(i[,1], i[,2], paste, sep=",")))
        
        ## select hap.p compatible with RP
        
        rows.comp <- lapply(hap.r.p, function(i) {
            tmp <- hap.p.no.gt2(s.M, M=M.all, i)
            tmp <- tmp[tmp>0]
            return(tmp)
        })
        

        ## get g.rsnp for haps compatible with RP
        h.rsnp <- lapply(rows.comp, function(i) sapply(names(i), function(j) as.numeric(sapply(strsplit(j, ","), function(k) substr(k, nchar(k), nchar(k))))))
        
        g.rsnp <- lapply(h.rsnp, function(i) {
            tmp <- colSums(i)
            inv <- i[1,tmp==1]==1 ## get hap 1 when rsnp is het
            tmp[names(tmp) %in% names(inv)[inv]] <- -1 ## if hap1 is 1, genotype is -1 (hap rsnp 1|0)
            return(tmp)
            })
        
        ## get p(H|G fsnps) for each hap pair in rows.comp
        pH.G <- lapply(rows.comp, function(i) apply(M.all.cond[i,,drop=F], 1, sum))

        ## prepare lists to return

        n.mat2 <- Reduce(cbind, lapply(1:ncol(n.mat), function(i) matrix(rep(n.mat[,i], length(g.rsnp[[i]])),ncol=length(g.rsnp[[i]]))))
        n.list <- lapply(1:nrow(n.mat), function(i) n.mat2[i,])
        g.list <- lapply(1:nrow(n.mat), function(i) unname(Reduce(c,g.rsnp)))
        pH.list <- lapply(1:nrow(n.mat), function(i) unname(Reduce(c,pH.G)))
        names(n.list) <- names(g.list) <- names(pH.list) <- rownames(n.mat)

        ## append

        n.all <- c(n.all, n.list)
        g.all <- c(g.all, g.list)
        pH.all <- c(pH.all, pH.list)
    }
    
    NB <- list(counts=y,p.g=p.g.comp)

    ## need to sort n.all, pH.all and g.all as in stan.f$m
    m=stan.f$m
    ASE <- list(m=m,n=n.all[names(m)], p=pH.all[names(m)], g=g.all[names(m)])

    return(list(NB=NB, ase=ASE))
    
    
}

           


#' sub-function for formatting input for trecase stan no GT, deals with genotype of rsnp and p(Grsnp|G fsnp) & p(H|G) to make stan input for negative binomial
#'
#' This function allows you to format input for stan Nb per individual: adapted from stan.rsnp.noGT.eff
#' @param counts.g DT of total gene counts, colnames sample names
#' @param M  matrix of haplotypes from reference panel for fsnps and rsnp, rows=hap.pairs, cols GT
#' @param M.cond  matrix haplotypes from reference panel, rows as M, cols GT, colsums=1
#' @param NB data table with summarize information for neg binom side, output from fsnp.prep or stan.fsnp.noGT
#' @keywords stan input rsnp
#' @export
#' 
#' @return list of stan input: list of 1) counts: data table with total counts  and 2) list p.g with p(genotype) rsnp per sample.
#' 
#' stan.nb()
#' 

stan.nb <- function(counts.g, M, M.cond, NB){

    ## make data table with no 0 entries in M to get genotype of rsnp for a given hap.pair
    gt.rf <- data.table(summary(M))
    setkey(gt.rf,i)    

    ## for each stan.f$NB$g.fsnps get p(gf+rsnp| gfsnp)
    u <- as.character(unique(NB$g.fsnps))
    p.g <- lapply(u,  function(i) {
        w <-  which(M.cond[,i]!=0)
        tmp <- M.cond[w,i]
        ## for each hap pair in names(tmp) get their corresponding genotype (colname with no-0 entry)
        gt.rf2 <- colnames(M)[gt.rf[i %in% w,j]]
        
        names(tmp) <- substr(gt.rf2,nchar(gt.rf2[1]), nchar(gt.rf2[1])) ## last character from gt.rf is geno rsnp
        tmp <- tapply(tmp, names(tmp), sum) ## add elements in vector with same genotype
        return(tmp)
    })
    names(p.g) <-u


    ## assign p.g to each sample,

    if(nrow(NB) == 1){
        
        samp.id <- lapply(1:nrow(NB), function(i) sort(unlist(strsplit(NB$ind[i], ","))))
    } else {
        samp.id  <- NB$ind 
    }
        
   
    p.g.comp <- lapply(1:length(samp.id), function(i) rep(list(p.g[[i]]), length(samp.id[[i]])) )

    p.g.comp <- Reduce(c,p.g.comp)  ## remove nested list

    ##  total counts for relevant samples, in same order as p.g.comp
   
    y <- counts.g[,unlist(samp.id),with=F]


    nb <- list(counts=y,p.g=p.g.comp)

    return(nb)
}
        


#' sub-function for formatting input for trecase stan no GT, deals with genotype of rsnp and p(Grsnp|G fsnp) & p(H|G) to make stan input for ASE
#'
#' This function allows you to format input for stan ASE per individual: adapted from stan.rsnp.noGT.eff
#' @param M  matrix of haplotypes from reference panel for fsnps and rsnp, rows=hap.pairs, cols GT
#' @param M.cond  matrix haplotypes from reference panel, rows as M, cols GT, colsums=1
#' @param s.M named vector with names rownames(M) and value swapped rownames(M), to identify swaps, done in stan.rsnp.noGT.eff 
#' @param n.mat matrix with hap for fsnps and ase counts, output from fsnp.prep or stan.fsnp.noGT
#' @param ai0.mat matrix with estimates for allelic imbnalance for reference bias correction, defaults to NULL
#'  @param vai0.mat matrix with variance estimates for allelic imbnalance for reference bias correction, defaults to NULL
#' @keywords stan input rsnp
#' @export
#' 
#' @return list of stan input: list of 1) list with vectors of n.counts, 2) list of vectors of  p(H(f+r)|Gf) , 3) list of vectors of possible genotype rsnp in scale 0,1,-1,2, each vector one sample, 4)optional list with AI reference bias estimates
#' 
#' stan.ase()
#' 

stan.ase <- function(M, M.cond, s.M, n.mat, ai0.mat=NULL, vai0.mat=NULL) {
    
    ## get hap pairs for fsnps and potentially any rsnp
    hap.f.p <- lapply(colnames(n.mat), function(i) unlist(strsplit(i,",")))
    hap.r.p <- lapply(hap.f.p, function(i) sapply(i, function(j) paste0(j,0:1)))
    hap.r.p <- lapply(hap.r.p, function(i) as.vector(outer(i[,1], i[,2], paste, sep=",")))
    
    ## select hap.p compatible with RP
        
    rows.comp <- lapply(hap.r.p, function(i) {
        tmp <- hap.p.no.gt2(s.M, M=M, i)
        tmp <- tmp[tmp>0]
        return(tmp)
    })
        

    ## get g.rsnp for haps compatible with RP
    h.rsnp <- lapply(rows.comp, function(i) sapply(names(i), function(j) as.numeric(sapply(strsplit(j, ","), function(k) substr(k, nchar(k), nchar(k))))))
    
    g.rsnp <- lapply(h.rsnp, function(i) {
        tmp <- colSums(i)
        inv <- i[1,tmp==1]==1 ## get hap 1 when rsnp is het
        tmp[names(tmp) %in% names(inv)[inv]] <- -1 ## if hap1 is 1, genotype is -1 (hap rsnp 1|0)
        return(tmp)
    })
        
    ## get p(H|G fsnps) for each hap pair in rows.comp
    pH.G <- lapply(rows.comp, function(i) rowSums(M.cond[i,,drop=F]))

    ## prepare lists to return

    n.mat2 <- Reduce(cbind, lapply(1:ncol(n.mat), function(i) matrix(rep(n.mat[,i], length(g.rsnp[[i]])),ncol=length(g.rsnp[[i]]))))
    n.list <- lapply(1:nrow(n.mat), function(i) n.mat2[i,])
    g.list <- lapply(1:nrow(n.mat), function(i) unname(Reduce(c,g.rsnp)))
    pH.list <- lapply(1:nrow(n.mat), function(i) unname(Reduce(c,pH.G)))
    names(n.list) <- names(g.list) <- names(pH.list) <- rownames(n.mat)

    
    ASE <- list(n=n.list, p=pH.list, g=g.list)

    if(!is.null(ai0.mat)){        
        ASE$ai0 <- aux.stan.ase(ai0.mat, g.rsnp, n.mat)
        ASE$vari0 <- aux.stan.ase(vai0.mat, g.rsnp, n.mat)
        
    }
    

    return(list(ASE))
}

#' aux function to prepare inputs for allelic imbalance estimates, wrapper to avoid code repetition
#' @param mat matrix with allelic imbalance or variance allelic imbalance estimates
#' @param g.rsnp object created in stan.ase containing genotypes compatible with reference panel
#' @param n.mat matrix with hap for fsnps and ase counts, output from fsnp.prep or stan.fsnp.noGT
#' @keywords aux reference panel input stan 
#' @export
#' @return list of stan input: list with AI reference bias estimates or variance estimates per individual
#' 
#' aux.stan.ase()

aux.stan.ase <- function(mat, g.rsnp, n.mat){
    ai0.mat <- Reduce(cbind, lapply(1:ncol(mat), function(i) matrix(rep(mat[,i], length(g.rsnp[[i]])),ncol=length(g.rsnp[[i]]))))
    ai0.list <- lapply(1:nrow(ai0.mat), function(i) ai0.mat[i,])
    names(ai0.list)  <- rownames(n.mat)
    return(ai0.list)

}
