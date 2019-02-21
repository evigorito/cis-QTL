library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggrepel)

#' get summaries from stan runs into a data table 
#'
#' This function allows you to combine files that can be opened as data tables into 1 data table, useful for stan summaries, tags, etc
#' @param path path to dir with files to combine
#' @param pattern pattern of files to combine
#' @keywords combine files into data table
#' @export
#' @return data table of combined files
#' comb.files()

comb.files <- function(path,pattern) {
    tmp <- list.files(path,pattern, full.names=TRUE)
    tmp <- rbindlist(lapply(tmp,fread))
    return(tmp)
}




#' combine summaries from stan models allowing different tags in 1 model
#'
#' This function allows you to combine summaries from diffrent stan models allowing for different tags in 1 model
#' @param x data table with summary for model 1, same tags as model 2
#' @param y data table with summary for model 2, same tags as model 1
#' @param z data table with summary for model 3, soem share some different tags
#' @param s suffixes for model 1,2,3, defaults to lm, gt, ngt
#' @param tags.m2 tags run in models 1 and 2: data table with Gene_id, tag, SNP
#' @param tags.m3 tags run in model 3: data table with Gene_id, tag, SNP
#' @keywords combine files into data table
#' @export
#' @return data table of combined files
#' comb.sum()

comb.sum <- function(x,y,z,s=c(".lm",".gt", ".ngt"),tags.m2,tags.m3) {

    ## get common genes in all summaries

    genes <- Reduce(intersect, list(unique(x$Gene_id), unique(y$Gene_id), unique(z$Gene_id)))

    x <- x[Gene_id %in% genes,]
    y <- y[Gene_id %in% genes,]
    z <- z[Gene_id %in% genes,]

    ## rename to ease merging and column id
    names(x)[!names(x) %in% c("Gene_id", "tag") ] <- paste0(names(x)[!names(x) %in% c("Gene_id", "tag") ], s[1])
    names(y)[!names(y) %in% c("Gene_id", "tag") ] <- paste0(names(y)[!names(y) %in% c("Gene_id", "tag") ], s[2])
    names(z)[!names(z) %in% c("Gene_id", "tag") ] <- paste0(names(z)[!names(z) %in% c("Gene_id", "tag") ], s[3])
    

    ## same tag
    x.y <- merge(x,y,by=c("Gene_id", "tag"), sort=F, suffixes=s[1:2])
    x.y.z <- merge(x.y,z, by=c("Gene_id", "tag"), sort=F, suffixes=s[2:3], all=FALSE)

    ## look for matching snp when tag is not shared.
    all <- merge(x.y,z,  by=c("Gene_id", "tag"), sort=F, suffixes=s[2:3], all=TRUE)
    ## get gene/tag pairs, first check if names are log2( or log2_
    n <- grep("se_mean.gt", names(all), value=T)
    if(length(grep("log2\\(", n))) {
        no.match.z <- all[is.na(`log2(aFC)_mean.ngt`) , .(Gene_id,tag)] ## tags run in gt/lm
        no.match.y <- all[is.na(`log2(aFC)_mean.gt`) , .(Gene_id,tag)] ## tags run in ngt
    } else {
        no.match.z <- all[is.na(log2_aFC_mean.ngt) , .(Gene_id,tag)] ## tags run in gt/lm
        no.match.y <- all[is.na(log2_aFC_mean.gt) , .(Gene_id,tag)] ## tags run in ngt

    }
    
    ## look for alternative snp in model z and y, respectively
    alt.z <- rbindlist(tag.comp(DT1=no.match.z, DT2=tags.m2, DT3=tags.m3, s=s[2]))
    alt.y <- rbindlist(tag.comp(DT1=no.match.y,DT2=tags.m3,DT3=tags.m2, s=s[3]))

    alt.zz <- merge(alt.z, z, by=c("Gene_id","tag"))
    ##setnames(alt.zz, names(alt.zz)[c(4:16,18:ncol(alt.zz))], paste0(names(alt.zz)[c(4:16,18:ncol(alt.zz))],s[3]))
    setnames(alt.zz,"tag",paste0("tag",s[3]))
    alt.z.all=merge(x.y, alt.zz, by.x=c("Gene_id","tag"), by.y=c("Gene_id",paste0("tag",s[2])))
    setnames(alt.z.all,"tag",paste0("tag", s[2]))

    alt.y <- merge(x.y, alt.y, by=c("Gene_id","tag"))
    setnames(alt.y,"tag",paste0("tag", s[2]))
    alt.y.all <- merge(alt.y, z,  by.x=c("Gene_id",paste0("tag", s[3])), by.y=c("Gene_id","tag"))#,suffixes=c("", s[3]))


    setcolorder(alt.z.all,names(alt.y.all))

    comb1 <- rbind(alt.y.all,alt.z.all)
    ##setnames(comb1, names(comb1)[c(11:23,25:26)], paste0(names(comb1)[c(11:23,25:26)], s[2]))

    ## rbind comb with x.y.z, lm.gt.no:make tag.gt same as tag and relabel tag as tag.ngt

    x.y.z[,tag.gt:=tag][,tag.ngt:=tag][,tag:=NULL]

    ## order columns as in comb1 to rbind
    setcolorder(x.y.z, names(comb1))

    comb <- rbind(comb1, x.y.z)
    ## add column as whether same tag was used in gt/lm vs nogt

    comb[get(paste0("tag",s[2]))==get(paste0("tag",s[3])), same.tag:="yes"][get(paste0("tag",s[2]))!=get(paste0("tag",s[3])), same.tag:="no"]

    return(comb) 
    
}

#' combine summaries from variable number of stan models allowing different tags in last  model
#'
#' This function allows you to combine summaries from many stan models allowing for different tags in 1 model
#' @param summ list with data tables with summary for each model, models with same tags except the last
#' @param prefix character vector with the prefixes of the tag and summary columns that identify the two key models 
#' @param s suffixes for models when merging to identify them
#' @param tags list with tags run in models 1-(n-1) and model n: each element a data table with Gene_id, tag, SNP
#' @keywords combine summaries tags
#' @export
#' @return data table of combined summaries
#' comb.sum.v()

comb.sum.v <- function(summ, prefix=c(".ngt", ""), s,tags) {

    ## get common genes in all summaries and tags

    genes <- Reduce(intersect, lapply(summ, function(i) unique(i$Gene_id)))
    summ <- lapply(summ, function(i) i[Gene_id %in% genes,])
    tags <- lapply(tags, function(i) {
        i[Gene_id %in% genes,]
        setkey(i, Gene_id, tags)
        return(i)
        })


   
    ## same tag

    tag.name <- paste0("tag",prefix)

    same.tag <- Reduce((function() { counter = 1
        function(x,y) {         
            d=merge(x,y,all=FALSE,by.x =c("Gene_id", tag.name[counter]), by.y =c("Gene_id", tag.name[counter+1]), sort=F, suffixes=s[counter:(counter+1)])
            counter <<- counter+2  ## super assignment operator to change counter in upper environment
            return(d)
            }})(), summ )
 
    ## look for matching snp when tag is not shared.
    ## merge many datatables
    all <- Reduce((function() { counter = 1
        function(x,y) {         
            d=merge(x,y,by.x =c("Gene_id", tag.name[counter]), by.y =c("Gene_id", tag.name[counter+1]), sort=F, suffixes=s[counter:(counter+1)],  all=TRUE)
            counter <<- counter+2
            return(d)
        }}) (), summ )


#### working in this function, not sure if worth it ###
    
    ## get gene/tag pairs
    no.match.z <- all[is.na(get(paste0("log2(aFC)_mean", prefix[1]))) , .(Gene_id,tag.ngt)] ## tags run rna
    no.match.y <- all[is.na(get(paste0("log2(aFC)_mean",prefix[2]))) , .(Gene_id,tag.ngt)] ## tags run ngt
    
    ## look for alternative snp in model z and y, respectively
    alt.z <- rbindlist(tag.comp(DT1=no.match.z, DT2=tags[[2]], DT3=tags[[1]], suf=prefix[1], s=s[2]))
    alt.y <- rbindlist(tag.comp(DT1=no.match.y,DT2=tags[[1]],DT3=tags[[2]], s=s[1]))

    alt.zz <- merge(alt.z, z, by=c("Gene_id","tag"))
    setnames(alt.zz, names(alt.zz)[c(4:16,18:ncol(alt.zz))], paste0(names(alt.zz)[c(4:16,18:ncol(alt.zz))],s[3]))
    setnames(alt.zz,"tag",paste0("tag",s[3]))
    alt.z.all=merge(x.y, alt.zz, by.x=c("Gene_id","tag"), by.y=c("Gene_id",paste0("tag",s[2])))
    setnames(alt.z.all,"tag",paste0("tag", s[2]))

    alt.y <- merge(x.y, alt.y, by=c("Gene_id","tag"))
    setnames(alt.y,"tag",paste0("tag", s[2]))
    alt.y.all <- merge(alt.y, z,  by.x=c("Gene_id",paste0("tag", s[3])), by.y=c("Gene_id","tag"),suffixes=c("", s[3]))


    setcolorder(alt.z.all,names(alt.y.all))

    comb1 <- rbind(alt.y.all,alt.z.all)
    setnames(comb1, names(comb1)[c(11:23,25:26)], paste0(names(comb1)[c(11:23,25:26)], s[2]))

    ## rbind comb with x.y.z, lm.gt.no:make tag.gt same as tag and relabel tag as tag.ngt

    x.y.z[,tag.gt:=tag][,tag.ngt:=tag][,tag:=NULL]

    ## order columns as in comb1 to rbind
    setcolorder(x.y.z, names(comb1))

    comb <- rbind(comb1, x.y.z)
    ## add column as whether same tag was used in gt/lm vs nogt

    comb[get(paste0("tag",s[2]))==get(paste0("tag",s[3])), same.tag:="yes"][get(paste0("tag",s[2]))!=get(paste0("tag",s[3])), same.tag:="no"]

    return(comb) 
    
}



            
    

#' get given a tag snp in one set it looks at if at least one of the tagged snps is part of a group in a different set and if so gives the gene, tag in set 1 and tag in set2. 
#'
#' This function allows you to compare outputs for models with tags
#' @param DT1 data table with Gene_id and tag in set1
#' @param DT2 data table with Gene_id,  tags and SNP in set1, extra cols allowed
#' @param DT3 data table with Gene_id tag and SNP for set2, extra cols allowed
#' @param suf suffix for tag col in DT1, defaults to NULL
#' @param s suffix for tag dataset1
#' @keywords stan output tag
#' @export
#' @return list of data tables per gene each with cols Gene_id tag in set1 tag set2
#' tag.comp()

tag.comp <- function(DT1,DT2,DT3,suf=NULL,s){

    setkey(DT2, Gene_id,tag,SNP)
    setkey(DT3, Gene_id,SNP)
    if(!is.null(suf)){
        tag.col <- paste0("tag",suf)
    } else{
        tag.col <- paste0("tag","")
    }
    
    

    alt <- mclapply(unique(DT1$Gene_id), function(i) {
        DT2.sub <- DT2[Gene_id==i,]
        DT3.sub <- DT3[Gene_id==i,]
        tgs=DT1[Gene_id==i, get(tag.col)]
        ## take the first SNP per gene that is tagged by tgs in gt and has a match in ngt.t.r
        snp <- lapply(tgs, function(j) DT2.sub[tag==j & SNP %in% DT3.sub$SNP, SNP][1])
        names(snp) <- tgs
        snp <- unlist(snp[!sapply(snp,is.na)])
        if(is.null(snp)){
            return(NULL)
        } else {
            
            tmp <- DT3.sub[SNP %in% snp,.(Gene_id,tag,SNP)]
            ## order tmp as snp
            tmp <- tmp[order(match(SNP,snp))]
            ## add tag
            tmp[,paste0("tag",s):=names(snp)]
            ##remove SNP
            tmp[,SNP:=NULL]
            return(tmp)
        } } )
        names(alt) <- unique(DT1$Gene_id)
        ## remove nulls
        alt <- alt[!sapply(alt,is.null)] 
        return(alt)
}


#' Correct direction of effect when runnning a SNP with a tag, based on EAF .
#'
#' This function allows you to compare account for deffinition of REF allele between a tag SNP and the SNP tagged in the direction of effects
#' @param DT data table with 2 tag cols plus stan output
#' @param le.file full path to legend file for the corresponding chr to look-up for eaf
#' @param s vector with suffixes for tag column names
#' @keywords stan output tag
#' @export
#' @return data table with extra cols with eaf for tags and corrected direction of effect size
#' eaf.aj()

eaf.aj <- function(DT,le.file,s){
    eaf.dif.tags <- apply(DT[same.tag=="no", paste0("tag",s), with=FALSE], 2, snp.eaf, file1='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz')

    eaf.dif.tags <- Reduce(cbind,eaf.dif.tags)
    names(eaf.dif.tags) <- unlist(lapply(s, function(i) paste0(c("tag","eaf"),i)))

    eaf.dif.tags[, op.dir:="no"][(eaf.gt<0.5 & eaf.ngt>0.5) | (eaf.gt>0.5 & eaf.ngt<0.5), op.dir:="yes"]
    ## remove duplicates

    setkeyv(eaf.dif.tags, names(eaf.dif.tags))
    eaf.dif.tags <- unique(eaf.dif.tags)

## some tag.ngt/tag.gt pairs have opposite allele freq. Include this info to assess direction of effects.

    comb <- merge(DT, eaf.dif.tags, by=paste0("tag",s), all.x=T)

    ## invert effect of "log2(aFC)_mean.ngt" when op.dir "yes"
     n <- grep("se_mean.gt", names(all), value=T)
    if(length(grep("log2\\(", n))) {
        comb[ op.dir=="yes", `log2(aFC)_mean.ngt`:=-`log2(aFC)_mean.ngt`]
    } else {
       comb[ op.dir=="yes", log2_aFC_mean.ngt:=-log2_aFC_mean.ngt]

    }

    

    return(comb)
}

#' Compare ngt stan output to lm and gt.
#'
#' Wrapping function to compare outputs when they are run with different tags, find a common tag and corrects effect size if eaf are coded in opposite directions
#' @param path.s to files with stan summaries for nogt model
#' @param pattern.s of summary files to read
#' @param path.t path to files with tag files for ngt model, defaults to NULL if same as path.s
#' @param pattern.t pattern for tag files
#' @param lm.sum data table with summary of lm
#' @param gt.sum data table with summary of gt
#' @param s vector with suffixes for tag column names for lm, gt and ngt
#' @param tags.m2 data table with tags run in model2 (gt)
#' @param le.file full path to legend file for the corresponding chr to look-up for eaf, use null if info already in output
#' @param 
#' @keywords stan output tag
#' @export
#' @return data table with extra cols with eaf for tags and corrected direction of effect size
#' comp.ngt()

comp.ngt <- function(path.s,pattern.s, path.t=NULL , pattern.t, lm.sum, gt.sum, s, tags.m2, le.file=NULL) {

    sum <- comb.files(path.s,pattern.s)
    
    if(is.null(path.t)){
        tags <- comb.files(path=path.s,pattern.t)
    }
    if(!is.null(path.t)) {
        tags <- comb.files(path.t,pattern.t)
    }

    
    ##tags run in model
    t.r <- merge(tags, sum[,.(Gene_id,tag)], by= c("Gene_id","tag"))

    ## combine with lm.g and gt.g
    tmp <- comb.sum(x=lm.sum,y=gt.sum,z= sum, s,tags.m2 , tags.m3= t.r)
    ## add EAF
    if(!is.null(le.file)){
        tmp <- eaf.aj(tmp,le.file, s[2:3])
    } else { ## correct direction of estimate
        tmp[, op.dir:='no'][(tag.EAF.gt<0.5 & tag.EAF.ngt>0.5) | (tag.EAF.gt>0.5 & tag.EAF.ngt<0.5), op.dir:="yes"]
        col.ngt <- grep("log2.*[^se]_mean.ngt", names(tmp), value=T)
        tmp[op.dir=='yes', (col.ngt) := -get(col.ngt) ]
    }
    
    
    return(unique(tmp))
}

#' Compare stan output per gene
#'
#' Compare number of sig and total SNPs called by gene between models
#' @param DT data table  with stan outputs, uses column `log2(aFC)_null` for each model
#' @keywords stan output compare SNPS by gene
#' @export
#' @return data table number of Significant and total SNPs per gene per model
#' snp.gene()

snp.gene <- function(DT){
    n <- grep("se_mean.gt", names(all), value=T)
    if(length(grep("log2\\(", n))) {
        cols <- grep("log2\\(aFC\\)_null", names(DT), value=T)
        s <- gsub("log2\\(aFC\\)_null", "", cols)
    } else {
        cols <- grep("log2_aFC_null", names(DT), value=T)
        s <- gsub("log2_aFC_null", "", cols)

    }
    
    
    tmp <- DT[,.N, c("Gene_id",cols )]  
    setnames(tmp,cols, paste0("null",s))
    ## convert to wide
    tmp <- t(dcast(tmp, paste(paste(paste0("null",s), collapse="+"), "~ Gene_id"),))
    ## format
    tmp <- tmp[-(1:2),]
    tmp[is.na(tmp)] <- 0
    rw.tot <- apply(tmp,1,function(i) sum(as.numeric(i)))
    dt <- data.table(apply(tmp,2,as.numeric))
    dt[,Gene_id:=rownames(tmp)]
    dt[,N:=rw.tot]
    setcolorder(dt,c("Gene_id", paste0("V", 1:4),"N"))
    setnames(dt,paste0("V", 1:4), c("sig.both","sig.ngt.null.rna","null.ngt.sig.rna","null.both"))
    setkey(dt,N)
    return(dt)
}

#' add Signif column to btrecase output comparing 2 conditions
#'
#' Codes the column with Both, None, option1, option2
#' @param dt data table with btrecase output
#' @param x1 name of cols with  null for condition1
#' @param x2 name of cols with  null for condition2
#' @param col identifier for first and second column being significant
#' @keywords Btrecase comparison significant  
#' @export
#' @return data table inout with new column
#' add.signif()

add.signif <- function(dt, x1, x2, col){

     null.cols <- c(x1, x2)
    dt[, Signif:="None"][get(null.cols[1])=="no" & get(null.cols[2])=="no", Signif:="Both"][get(null.cols[1])=="no" & get(null.cols[2])=="yes", Signif:=col[1]][get(null.cols[1])=="yes" & get(null.cols[2])=="no", Signif:=col[2]]
    dt[,Signif:=factor(Signif, levels=c("None", col, "Both"))]
    return(dt)

}




#' plots a parameter ditribution from Btrecase output
#'
#' This function allows you to  compare parameter estimates from 2 conditions
#' @param dt data table with parameter mean, CI min and max
#' @param x1 names of cols with mean, CI min and max and null for condition1
#' @param x2 names of cols with mean, CI min and max and null for condition2
#' @param title optional title
#' @param rx vector with range for x axis
#' @param ry vector with range for y axis
#' @param xl character with xlab
#' @param yl character with ylab
#' @param col when usimg color for symbol based on sig or not sig, provide identifier for first and second column being significant
#' @param title plot
#' @param axis.title size for axis.title defaults to 30
#' @param axis.text size for axis.text.x and axis.text.y defaults to 28
#' @param legend.title size for legend title, defaults to 22
#' @param legend.text size for legend.text, defaults to 18
#' @param legend.symbol size for legend symbols, defaults to 10
#' @param point size for geom_points, defaults to 6
#' @keywords Btrecase summary plot  
#' @export
#' @return ggplot object
#' btrecase.plot()

btrecase.plot <- function(dt, x1, x2, t=NULL, e=NULL, rx=NULL,ry=NULL,xl, yl, col=NULL, title=NULL,
                          axis.title=30, axis.text=28, legend.title=22, legend.text=18, legend.symbol=10, point.size=6){   
    old.n <- function(x,pat){
        tmp <- grep(pat,x,value=T)
        return(tmp)
    }
    sub.grep <- function(x, pat, sub){
        tmp <- old.n(x,pat)
        tmp2 <- unname(sapply(tmp, function(i) gsub(pat, sub, i)))
        return(tmp2)
    }
    x <- copy(dt)
    if(nrow(x)> 5000) { ## sample 5000
        x <- x[sample(nrow(x), 5000),]
    }
    
    new.n <-  sub.grep(names(x), "%", "")
    setnames(x, old.n(names(x), "%"), new.n)
    x1 <- sapply(x1, gsub, pattern="%",replacement= "", USE.NAMES=F)
    x2 <- sapply(x2, gsub, pattern="%",replacement= "", USE.NAMES=F)

    if(any(names(x)=="op.dir")){
        ## need to make CI quantiles compatible with estimate
        q.n <- new.n[new.n %in% x2]
        x[, paste0(q.n,".comp") := lapply(q.n, function(i) x[[i]])]
        x[op.dir=='yes', paste0(q.n,".comp") := lapply(rev(q.n), function(i) -unlist(x[op.dir=='yes',i,with=F]))]
        x[, (q.n) := lapply(paste0(q.n,".comp"), function(i) x[[i]])]
    }
       
    p <- ggplot(x, aes(get(x1[1]), get(x2[1]), label=Gene_id)) +      
        geom_errorbarh(aes_string(xmin=x1[2], xmax=x1[3]), linetype=2, colour='grey' ,size=0.1)  +
        geom_errorbar(aes(ymin=get(x2[2]), ymax=get(x2[3])),linetype=2, colour='grey', width=0.1 ,size=0.1)  +
        theme_bw() +
        xlab(xl) + ylab(yl) +
        theme(axis.title = element_text(size=axis.title), axis.text.x = element_text(colour="black", size = axis.text),axis.text.y = element_text(colour="black", size = axis.text), legend.title = element_text(size=legend.title), 
              legend.text = element_text(size=legend.text)) +
        geom_rug(col=rgb(.7,0,0,alpha=.2)) +
        geom_vline(xintercept=0,color = "blue") +
        geom_hline(yintercept=0,color = "blue") +
        geom_abline(intercept=0,slope=1, linetype=2, size=0.5) +
        guides(colour=guide_legend(override.aes = list(size = legend.symbol)))

    if(!is.null(col)){
        null.cols <- c(x1[4], x2[4])
        
        x[, Signif:="None"][get(null.cols[1])=="no" & get(null.cols[2])=="no", Signif:="Both"][get(null.cols[1])=="no" & get(null.cols[2])=="yes", Signif:=col[1]][get(null.cols[1])=="yes" & get(null.cols[2])=="no", Signif:=col[2]]
        x[,Signif:=factor(Signif, levels=c("None", col, "Both"))]
        
        p <- p + geom_point(aes(colour=Signif), shape=1, size=point.size) ##+
        ##geom_text_repel(aes(label=ifelse(get(x1[1]) %in% x[Gene_id %in% fpg, get(x1[1])], Gene_id, "")))
    } else {   
        p <- p + geom_point(shape=1, size=point.size)
    }

    
        
    if(!is.null(title)){
        p <- p + ggtitle(title) +
            theme(plot.title=element_text(size=32,hjust = 0.5))
        
    }
    if(!is.null(rx)){
        p <- p + scale_x_continuous(limits = rx)
    }
    if(!is.null(ry)){
        p <- p + scale_y_continuous(limits = ry)
    }
    
    return(p)
}

#' wrap functions to plot a parameter ditribution from stan summaries comparing 2 conditions
#'
#' This function allows you to  plot a parameter ditribution (output from stan.to.plot)
#' @param x list of stan summaries
#' @param y list of stan summaries
#' @param par parameter to extract from both list of summaries
#' @param sx suffix for names of x to distinguish from  y
#' @param sy suffix for names of y to distinguish from  x
#' @param t optional title
#' @param e simulated effect size to draw vertical and horizontal lines
#' @param rx vector with range for x axis
#' @param ry vector with range for x axis
#' @param xl character xlab
#' @param yl character yl
#' @keywords stan summary plot  
#' @export
#' @return ggplot object
#' wrap.plot()

wrap.plot <- function(x,y,param="bj",sx="_fix",sy="_prob",t,e=NULL, rx=NULL,ry=NULL,xl="Fixed haplotyes", yl="Uncertain haplotypes"){
    DTx <- stan.to.plot(x=x,y=param)
    DTy <- stan.to.plot(x=y, y=param)
    names(DTx) <- paste0(names(DTx), sx)
    names(DTy) <- paste0(names(DTy),sy)
    DT.plot <- stan.plot(x=DTx,y=DTy,t=t,e=e, rx=rx, ry=ry,xl,yl)
    return(DT.plot)
}

#' Plot to comapre gene-snps association between gt and ngt for common genes and snps
#'
#' This function allows you to  plot gene-snps associations between gt and nogt from snps matched by tags in both conditions
#' @param x data table with stan output with both gt and nogt, output of comp.ngt or var.e
#' @param y  numeric value for r2 to filter by
#' @param gene character vector with gene ID
#' @keywords stan plot gene-snps compare
#' @export
#' @return ggplot object
#' gene.plot()

gene.plot <- function(x,y=NULL,gene){
    ## select entries for gene
    dt <- x[Gene_id==gene,]
    ## create col for position, use tag.ngt because it has more entries
    ##dt[, pos.gt:=as.numeric(gsub(":.*","", tag.gt))]
    dt[, Position:=as.numeric(gsub(":.*","", tag.ngt))]
    ## order `log2(aFC)_null.gt` to make consistent coloured plots when only "yes" is present
    dt[, `log2(aFC)_null.gt`:=factor(`log2(aFC)_null.gt`, levels=c("yes", "no"))]
    dt[, `log2(aFC)_null.ngt`:=factor(`log2(aFC)_null.ngt`, levels=c("yes", "no"))]
    ## plot abs(`log2(aFC)_mean.xx), modify cred intervals when necessary

    dt[,abs.ci.low:=ifelse(`log2(aFC)_mean.gt`> 0, `log2(aFC)_2.5%.gt`, -`log2(aFC)_2.5%.gt`)]
    dt[,abs.ci.h:=ifelse(`log2(aFC)_mean.gt`> 0, `log2(aFC)_97.5%.gt`, -`log2(aFC)_97.5%.gt`)]
    dt[,abs.ci.ngt.low:=ifelse(`log2(aFC)_mean.ngt`>0, `log2(aFC)_2.5%.ngt`, -`log2(aFC)_2.5%.ngt`)]
    dt[,abs.ci.ngt.h:=ifelse(`log2(aFC)_mean.ngt`>0, `log2(aFC)_97.5%.ngt`, -`log2(aFC)_97.5%.ngt`)]

     ## make plots with same x axis and same y axis
    min.pos <- min(dt$Position)
    max.pos <- max(dt$Position)

    min.y.axis <- round(min(c(dt$abs.ci.low, dt$abs.ci.ngt.low)) - 0.1 , 2)
    max.y.axis <- round(max(c(dt$abs.ci.h, dt$abs.ci.ngt.h)) + 0.1, 2)

    if(!is.null(y)){
        dt <- dt[r2>=y,]
    }
    

    p <- ggplot(dt, aes(x=Position, abs(`log2(aFC)_mean.gt`))) + geom_point(aes(colour=`log2(aFC)_null.gt`)) +  geom_hline(yintercept=0) + geom_errorbar(data=dt, aes(ymin=abs.ci.low, ymax=abs.ci.h), color="grey",linetype="dashed",width=0.1) + ggtitle(gene) + ylab(label="|log2(aFC)|, known GT") + xlim(min.pos-10,max.pos+10) + ylim(min.y.axis, max.y.axis)

    p2 <- ggplot(dt, aes(x=Position, abs(`log2(aFC)_mean.ngt`))) + geom_point(aes(colour=`log2(aFC)_null.ngt`)) +  geom_hline(yintercept=0) + geom_errorbar(data=dt, aes(ymin=abs.ci.ngt.low, ymax=abs.ci.ngt.h), color="grey",linetype="dashed",width=0.1)+ ylab(label="|log2(aFC)|, unknown GT") + xlim(min.pos-10,max.pos+10) + ylim(min.y.axis, max.y.axis)

    a <- plot_grid(p,p2, ncol=1)
    return(a)
}


#' Plot to comapre gene-snps association between gt and ngt from all snps tested in each condition
#'
#' This function allows you to  plot all tested associations for a given gene comparing gt and nogt 
#' @param x data table with stan output for gt
#' @param y data table with stan output for ngt, output of comp.ngt or var.e or other
#' @param gene character vector with gene ID
#' @param suf suffix to use when different to "ngt"
#' @keywords stan plot gene-snps compare
#' @export
#' @return ggplot object
#' gene.plot2()

gene.plot2 <- function(x,y,gene,suf=NULL){
    ## select entries for gene
    x <- x[Gene_id==gene,]
    y <- y[Gene_id==gene,]
    ## create col for position for gt and ngt
    x[, Position:=as.numeric(gsub(":.*","", tag))]
    tag.y <- grep("tag",names(y),value=T)
    ## use first one
    y[, Position:=as.numeric(gsub(":.*","", get(tag.y[1])))]
    ## order `log2(aFC)_null.gt` to make consistent coloured plots when only "yes" is present
    x[, `log2_aFC_null`:=factor(`log2_aFC_null`, levels=c("yes", "no"))]
    if(is.null(suf)){
        y[, `log2_aFC_null.ngt`:=factor(`log2_aFC_null.ngt`, levels=c("yes", "no"))]
        min.y.axis <- min(c(x$`log2_aFC_2.5%`, y$`log2_aFC_2.5%.ngt`))
        max.y.axis <- max(c(x$`log2_aFC_97.5%`, y$`log2_aFC_97.5%.ngt`))        
    } else {
         y[, `log2(aFC)_null.ngt`:=factor(paste0("log2(aFC)_null",suf), levels=c("yes", "no"))]
        min.y.axis <- min(c(x$`log2_aFC_2.5%`, y[[paste0("log2_aFC_2.5%",suf)]]))
        max.y.axis <- max(c(x$`log2_aFC_97.5%`, y[[paste0("log2_aFC_97.5%",suf)]]))
    }
    
    ## make plots with same x axis and same y axis
    min.pos <- min(c(x$Position, y$Position))
    max.pos <- max(c(x$Position, y$Position))


    p <- ggplot(x, aes(x=Position, `log2_aFC_mean`)) +
        geom_point(aes(colour=`log2_aFC_null`)) +
        geom_hline(yintercept=0) +
        geom_errorbar(data=x, aes(ymin=`log2_aFC_2.5%`, ymax=`log2_aFC_97.5%`), color="grey",linetype="dashed",width=0.1) +
        ggtitle(gene) +
        ylab(label="log2(aFC), known GT") +
        xlim(min.pos-10,max.pos+10) + ylim(min.y.axis-0.05, max.y.axis+0.05)

    p2 <- ggplot(y, aes(x=Position, get(paste0("log2_aFC_mean",suf)))) + geom_point(aes(colour=get(paste0("log2_aFC_null",suf)))) +
        geom_hline(yintercept=0) +
        geom_errorbar(data=y, aes(ymin=get(paste0("log2_aFC_2.5%",suf)), ymax=get(paste0("log2_aFC_97.5%",suf))), color="grey",linetype="dashed",width=0.1) +
        ylab(label="log2(aFC), imputed GT")+
        xlim(min.pos-10,max.pos+10) +
        ylim(min.y.axis-0.05, max.y.axis+0.05)

    a <- plot_grid(p,p2, ncol=1)
    return(a)
}

#' Plot to comapre gene-snps association between gt, ngt and rna from all snps tested in each condition
#'
#' This function allows you to  plot all tested associations for a given gene comparing gt,nogt and rna gt of exonic snps
#' @param x data table with stan output for gt
#' @param ci.x  numeric with probability of ci containing true value for gt, defaults to 0.95 (2.5 and 97.5)
#' @param y data table with stan output for ngt and rna, output of comp.ngt or var.e merged with rna
#' @param ci.y  numeric with probability of ci containing true value for ngt, defaults to 0.95 (2.5 and 97.5)
#' @param ci.z  numeric with probability of ci containing true value for rna, defaults to 0.95 (2.5 and 97.5)
#' @param z data table with stan output for all combined to indicate comparable snps
#' @param info.s optional to add info score on second plot, give name of column with info
#' @param fisher optional to add lowest pvalue for fisher test of p(het) between RNA and DNA, give name of column with fisher test
#' @param gene character vector with gene ID
#' @keywords stan plot gene-snps compare
#' @export
#' @return ggplot object
#' gene.plot3()

gene.plot3 <- function(x,ci.x=NULL, y,ci.y=NULL, ci.z=NULL,z,gene,info.s=NULL, fisher=NULL){
    ## select entries for gene
    x <- x[Gene_id==gene,]
    y <- y[Gene_id==gene,]
    ## create col for position for gt and ngt
    x[, Position:=as.numeric(gsub(":.*","", tag))]
    ## get col names for "tag" in y DT
    tag.y <- grep("tag",names(y),value=T)
    ## use first one
    y[, Position:=as.numeric(gsub(":.*","", get(tag.y[1])))]
    ## order `log2(aFC)_null.gt` to make consistent coloured plots when only "yes" is present
    null <- grep("null", names(x),value=T)
    ## merge x with z to get common snps, null column for x from z
    x <- merge(x, z[,c('Gene_id', 'tag.gt', 'tag.ngt', paste0(null,'.gt')), with=F], by.x=c('Gene_id','tag'), by.y=c('Gene_id', 'tag.gt') , all.x=T)
    x[is.na(get(paste0(null, '.gt'))) , paste0(null, ".gt") := paste0(get(null), '.gt')]
    x[, paste0(null, ".gt"):=factor(get(paste0(null, ".gt")), levels=c("yes", "no", "yes.gt", "no.gt"))]
    x[,(null):= get(paste0(null, ".gt"))]

    x[, paste0(null, ".gt"):=NULL]

    ## for snps in z (comparable SNPs) change position to nogt tag to ease comparison in plots
    x[!is.na(tag.ngt), Position:=as.numeric(gsub(":.*","", tag.ngt))]
    

    ## get suffixes in y
    suf <- gsub('.*_null', "", grep('.*_null', names(y), value=T))
    ## get prefixes for null in y
    pre <- unique(mapply(function(x,y) {gsub(y,"",x)}, grep('.*_null', names(y), value=T),suf,USE.NAMES=F))

    ## merge y with z
    y <- merge(y,z[,c('Gene_id', 'tag.ngt', paste0(pre,c('.ngt','.rna'))), with=F], by.x=c('Gene_id','tag'), by.y=c('Gene_id', 'tag.ngt') , all.x=T, suffixes=c('',2))

    suf2 <- '.nogt'
    
    y[is.na(get(paste0(pre, c('.ngt','.rna2')))) , paste0(pre, c('.ngt','.rna2')) := lapply(1:2, function(i) paste0(get(paste0(pre, suf[i])), suf2))]

    suf3 <- c('.ngt','.rna2')
    
    y[,  paste0(pre, suf3):=lapply(1:2, function(i) factor(get(paste0(pre,suf3[i])) , levels=c("yes",  "no" , paste0(c('yes', 'no'), suf2)) ))]

    y[, (paste0(pre,suf)) := lapply(1:2, function(i) get(paste0(pre,suf3[i]))) ]
   
    ## make plots with same x axis and same y axis
    min.pos <- min(c(x$Position, y$Position))
    max.pos <- max(c(x$Position, y$Position))


    ## get prefixes for x and y
    pre.x <- gsub('_null', "", grep('.*_null', names(x), value=T))
    pre.y <- gsub('_null', "", pre)

    ## get ci limits
    p.x <- ifelse(is.null(ci.x), 0.95, ci.x)
    p.y <- ifelse(is.null(ci.y), 0.95, ci.y)
    p.z <- ifelse(is.null(ci.z), 0.95, ci.z)
    p.yz <- c(p.y,p.z)

    min.y.axis <- min(c(x[[paste0(pre.x,'_',(1-p.x)*50,'%')]], sapply(1:length(suf), function(i) y[[paste0(pre.y,  '_',(1-p.yz[i])*50,'%' , suf[i])]])))
    max.y.axis <- max(c(x[[paste0(pre.x,'_',(p.x+(1-p.x)/2)*100, '%')]], sapply(1:length(suf), function(i) y[[paste0(pre.y, '_',(p.yz[i] +(1-p.yz[i])/2)*100,'%' , suf[i])]])))

    colors <- c("red", "blue4", "lightpink3", "skyblue4")
    names(colors) <- c("yes", "no", "yes.gt", "no.gt")
    man.col <- colors[names(colors) %in% x[[null]]]

    
    
    p <- ggplot(x, aes(x=Position, get(paste0(pre.x,'_mean')), color=get(paste0(pre.x,'_null'))))  +  geom_hline(yintercept=0) + geom_errorbar(data=x, aes(ymin=get(paste0(pre.x,'_',(1-p.x)*50,'%')), ymax=get(paste0(pre.x,'_', (p.x+(1-p.x)/2)*100, '%'))), color="grey90",linetype="dashed",width=0.1) + ggtitle(gene)  + xlim(min.pos-10,max.pos+10) + ylim(min.y.axis-0.05, max.y.axis+0.05)+ labs(y="log2(aFC),GT", color=paste0(pre.x, "_null"))+geom_point(shape=1, size=3) + scale_color_manual(values=man.col)

    cond <- c("DNA pha-SNPs", "RNA pha-SNPs")
    names(colors) <- c("yes",  "no" , paste0(c('yes', 'no'), suf2))
    
    col.cond <- lapply(suf, function(i) colors[names(colors) %in% y[[paste0(pre.y, "_null", i)]]])
    p23 <- lapply(1:length(suf),  function(i) {
        p <- ggplot(y, aes(x=Position, get(paste0(pre.y,"_mean", suf[i])), color=get(paste0(pre.y,"_null", suf[i]))))  +  geom_hline(yintercept=0) +
            geom_errorbar(data=y, aes(ymin=get(paste0(pre.y, '_',(1-p.yz[i])*50,'%', suf[i])), ymax=get(paste0(pre.y,'_',(p.yz[i] +(1-p.yz[i])/2)*100,'%', suf[i]))), color="grey90",linetype="dashed",width=0.1) +
            xlim(min.pos-10,max.pos+10) +
            ylim(min.y.axis-0.05, max.y.axis+0.05) +
            labs(y=paste0("log2(aFC),",cond[i]), color=paste0(pre.y, "_null \n",cond[i])) + scale_color_manual(values=col.cond[[i]]) + geom_point(shape=1, size=3)
        if(i==1 & !is.null(info.s)){
                    labs <- paste("Median info  =" , round(mean(y[['info.ngt']]),1),'\u00B1', round(sd(y[['info.ngt']]),1))
                    p <- p + annotate(geom="text", label= labs, y = c(max.y.axis), x = min(x$Position)+ 0.01*( max(x$Position) - min(x$Position)), hjust=0, vjust=0)

        }
        if(i==2 & !is.null(fisher)){
                    labs <- paste("min p-value phasing-SNP  =" , formatC(unique(y[[fisher]]),digits=1 ))
                    p <- p + annotate(geom="text", label= labs, y = c(max.y.axis), x = min(x$Position)+ 0.01*( max(x$Position) - min(x$Position)), hjust=0, vjust=0)

        }
        
        
        return(p)
        })
    a <- plot_grid(plotlist= c(list(p),p23), ncol=1)
    return(a)
}

#' Plot to comapre gene-snps association between gt,and rna from all snps tested in each condition
#'
#' This function allows you to  plot all tested associations for a given gene comparing gt,nogt and rna gt of exonic snps
#' @param x data table with stan output for gt
#' @param ci.x  numeric with probability of ci containing true value for gt, defaults to 0.95 (2.5 and 97.5)
#' @param y data table with stan output for ngt and rna, output of comp.ngt or var.e merged with rna
#' @param ci.z  numeric with probability of ci containing true value for rna, defaults to 0.95 (2.5 and 97.5)
#' @param z data table with stan output for all combined to indicate comparable snps
#' @param info.s optional to add info score on second plot, give name of column with info
#' @param fisher optional to add lowest pvalue for fisher test of p(het) between RNA and DNA, give name of column with fisher test
#' @param gene character vector with gene ID
#' @keywords stan plot gene-snps compare
#' @export
#' @return ggplot object with gt and rna plots only
#' gene.plot2b()

gene.plot2b <- function(x,ci.x=NULL, y, ci.y=NULL, ci.z=NULL,z,gene,info.s=NULL, fisher=NULL){
    ## select entries for gene
    ## select entries for gene
    x <- x[Gene_id==gene,]
    y <- y[Gene_id==gene,]
    ## create col for position for gt and ngt
    x[, Position:=as.numeric(gsub(":.*","", tag))]
    ## get col names for "tag" in y DT
    tag.y <- grep("tag",names(y),value=T)
    ## use first one
    y[, Position:=as.numeric(gsub(":.*","", get(tag.y[1])))]
    ## order `log2(aFC)_null.gt` to make consistent coloured plots when only "yes" is present
    null <- grep("null", names(x),value=T)
    ## merge x with z to get common snps, null column for x from z
    x <- merge(x, z[,c('Gene_id', 'tag.gt', 'tag.ngt', paste0(null,'.gt')), with=F], by.x=c('Gene_id','tag'), by.y=c('Gene_id', 'tag.gt') , all.x=T)
    x[is.na(get(paste0(null, '.gt'))) , paste0(null, ".gt") := paste0(get(null), '.gt')]
    x[, paste0(null, ".gt"):=factor(get(paste0(null, ".gt")), levels=c("yes", "no", "yes.gt", "no.gt"))]
    x[,(null):= get(paste0(null, ".gt"))]

    x[, paste0(null, ".gt"):=NULL]

    ## for snps in z (comparable SNPs) change position to nogt tag to ease comparison in plots
    x[!is.na(tag.ngt), Position:=as.numeric(gsub(":.*","", tag.ngt))]
    

    ## get suffixes in y
    suf <- gsub('.*_null', "", grep('.*_null', names(y), value=T))
    ## get prefixes for null in y
    pre <- unique(mapply(function(x,y) {gsub(y,"",x)}, grep('.*_null', names(y), value=T),suf,USE.NAMES=F))

    ## merge y with z
    y <- merge(y,z[,c('Gene_id', 'tag.ngt', paste0(pre,c('.ngt','.rna'))), with=F], by.x=c('Gene_id','tag'), by.y=c('Gene_id', 'tag.ngt') , all.x=T, suffixes=c('',2))

    suf2 <- '.nogt'
    
    y[is.na(get(paste0(pre, c('.ngt','.rna2')))) , paste0(pre, c('.ngt','.rna2')) := lapply(1:2, function(i) paste0(get(paste0(pre, suf[i])), suf2))]

    suf3 <- c('.ngt','.rna2')
    
    y[,  paste0(pre, suf3):=lapply(1:2, function(i) factor(get(paste0(pre,suf3[i])) , levels=c("yes",  "no" , paste0(c('yes', 'no'), suf2)) ))]

    y[, (paste0(pre,suf)) := lapply(1:2, function(i) get(paste0(pre,suf3[i]))) ]
   
    ## make plots with same x axis and same y axis
    min.pos <- min(c(x$Position, y$Position))
    max.pos <- max(c(x$Position, y$Position))


    ## get prefixes for x and y
    pre.x <- gsub('_null', "", grep('.*_null', names(x), value=T))
    pre.y <- gsub('_null', "", pre)

    ## get ci limits
    p.x <- ifelse(is.null(ci.x), 0.95, ci.x)
    p.y <- ifelse(is.null(ci.y), 0.95, ci.y)
    p.z <- ifelse(is.null(ci.z), 0.95, ci.z)
    p.yz <- c(p.y,p.z)

    min.y.axis <- min(c(x[[paste0(pre.x,'_',(1-p.x)*50,'%')]], sapply(1:length(suf), function(i) y[[paste0(pre.y,  '_',(1-p.yz[i])*50,'%' , suf[i])]])))
    max.y.axis <- max(c(x[[paste0(pre.x,'_',(p.x+(1-p.x)/2)*100, '%')]], sapply(1:length(suf), function(i) y[[paste0(pre.y, '_',(p.yz[i] +(1-p.yz[i])/2)*100,'%' , suf[i])]])))

    colors <- c("red", "blue4", "lightpink3", "skyblue4")
    names(colors) <- c("yes", "no", "yes.gt", "no.gt")
    man.col <- colors[names(colors) %in% x[[null]]]

    
    
    p <- ggplot(x, aes(x=Position, get(paste0(pre.x,'_mean')), color=get(paste0(pre.x,'_null'))))  +  geom_hline(yintercept=0) + geom_errorbar(data=x, aes(ymin=get(paste0(pre.x,'_',(1-p.x)*50,'%')), ymax=get(paste0(pre.x,'_', (p.x+(1-p.x)/2)*100, '%'))), color="grey90",linetype="dashed",width=0.1) + ggtitle(gene)  + xlim(min.pos-10,max.pos+10) + ylim(min.y.axis-0.05, max.y.axis+0.05)+ labs(y="log2(aFC),GT", color=paste0(pre.x, "_null"))+geom_point(shape=1, size=3) + scale_color_manual(values=man.col)

    cond <- c("DNA pha-SNPs", "imputed")
    names(colors) <- c("yes",  "no" , paste0(c('yes', 'no'), suf2))
    
    col.cond <- lapply(suf, function(i) colors[names(colors) %in% y[[paste0(pre.y, "_null", i)]]])
    p23 <- lapply(1:length(suf),  function(i) {
        p <- ggplot(y, aes(x=Position, get(paste0(pre.y,"_mean", suf[i])), color=get(paste0(pre.y,"_null", suf[i]))))  +  geom_hline(yintercept=0) +
            geom_errorbar(data=y, aes(ymin=get(paste0(pre.y, '_',(1-p.yz[i])*50,'%', suf[i])), ymax=get(paste0(pre.y,'_',(p.yz[i] +(1-p.yz[i])/2)*100,'%', suf[i]))), color="grey90",linetype="dashed",width=0.1) +
            xlim(min.pos-10,max.pos+10) +
            ylim(min.y.axis-0.05, max.y.axis+0.05) +
            labs(y=paste0("log2(aFC),",cond[i]), color=paste0(pre.y, "_null \n",cond[i])) + scale_color_manual(values=col.cond[[i]]) + geom_point(shape=1, size=3)
        if(i==1 & !is.null(info.s)){
                    labs <- paste("Median info  =" , round(mean(y[['info.ngt']]),1),'\u00B1', round(sd(y[['info.ngt']]),1))
                    p <- p + annotate(geom="text", label= labs, y = c(max.y.axis), x = min(x$Position)+ 0.01*( max(x$Position) - min(x$Position)), hjust=0, vjust=0)

        }
        if(i==2 & !is.null(fisher)){
                    labs <- paste("min p-value phasing-SNP  =" , formatC(unique(y[[fisher]]),digits=1 ))
                    p <- p + annotate(geom="text", label= labs, y = c(max.y.axis), x = min(x$Position)+ 0.01*( max(x$Position) - min(x$Position)), hjust=0, vjust=0)

        }
             
        return(p)
    })
    
    a <- plot_grid(plotlist= c(list(p,p23[[2]])), ncol=1)
    return(a)
}



#' Plot to comapre gene-snps association between gt, ngt and rna for the intersenction of the snps tested in each condition
#'
#' This function allows you to  plot common associations for a given gene comparing gt,nogt and rna gt of exonic snps
#' @param DT data table with stan output for gt, ngt and rna
#' @param gene character vector with gene ID
#' @param info.s optional to add info score on second plot, give name of column with info
#' @param fisher optional to add lowest pvalue for fisher test of p(het) between RNA and DNA, give name of column with fisher test
#' @keywords stan plot gene-snps compare
#' @export
#' @return ggplot object
#' gene.plot4()

gene.plot4 <- function(DT,gene, info.s=NULL, fisher=NULL){
    
    ## select entries for gene
    x <- DT[Gene_id==gene,]
    ## create col for position 
    x[, Position:=as.numeric(gsub(":.*","", tag.gt))]
    ## remove lm cols
    rem <- grep("lm",names(x),value=T)
    x[, (rem):=NULL]
    ## get suffixes
    suf <- gsub('.*_null', "", grep('.*_null', names(x), value=T))
    ## get prefixes for null 
    pre <- unique(mapply(function(x,y) {gsub(y,"",x)}, grep('.*_null', names(x), value=T),suf,USE.NAMES=F))
    x[, paste0(pre, suf):=lapply(suf, function(i) factor(get(paste0(pre,i)) , levels=c("yes", "no")))]
    ## get prefixes more general
    pre.x <- gsub('null', "", pre)
    ## get null cols
    null <- grep("null", names(x),value=T)
    colors <- c("red", "blue4")
    names(colors) <- c("yes", "no")
 
    cond <- c("DNA GT", "DNA pha-SNPs", "RNA pha-SNPs")
    col.cond <- lapply(null, function(i) colors[names(colors) %in% x[[i]]])

    ## correct CI before plotting
    q.n <- lapply(suf[2:3], function(i) sapply(c("2.5%", "97.5%"), function(j) paste0(pre.x, j,i), USE.NAMES=F))
    x <- lapply(q.n, function(i) {
        x[, paste0(i,".comp") := lapply(i, function(j) x[[j]])]
        x[op.dir=='yes', paste0(i,".comp") := lapply(rev(i), function(j) -unlist(x[op.dir=='yes',j,with=F]))]
        x[, (i) := lapply(paste0(i,".comp"), function(j) x[[j]])]
    })[[1]] ## get first elements, they are equal as x was modified in place

    min.y.axis <- min(x[,paste0(pre.x,"2.5%",suf), with=F])
    max.y.axis <- max(x[,paste0(pre.x,"97.5%",suf), with=F])

    
    p <- lapply(1:length(suf),  function(i) {
            p <- ggplot(x, aes(x=Position, get(paste0(pre.x,"mean", suf[i])), color=get(null[i]))) +  geom_hline(yintercept=0) +
                geom_errorbar(data=x, aes(ymin=get(paste0(pre.x, "2.5%", suf[i])), ymax=get(paste0(pre.x,"97.5%", suf[i]))), color="grey",linetype="dashed", width=0 ,size=0.1) +
                ylim(min.y.axis-0.05, max.y.axis+0.05) +
                labs(y=paste0("log2(aFC),",cond[i]), color=paste0(pre.x, "null \n",cond[i])) + scale_color_manual(values=col.cond[[i]]) + geom_point()
            if(i==2 & !is.null(info.s)){
                if(nrow(x) <=30 ) {
                ## remove duplicate values for info in consecutive positions
                x <- setkey(x, Position)
                x[,info2:=as.character(round(get(info.s),1))]
                rle.inf <- rle(x$info2)
                w <- which(rle.inf$length > 1) ## repeated consecutive elements
                l <- rle.inf$length
                for(j in w){
                    if(j==1) { ## first element repeated
                        x$info2[2: l[1] ]   <- ""
                    } else {
                        x$info2[(sum(l[1:(j-1)]) + 2) : sum(l[1:j]) ]   <- ""
                        }                 
                }
                
                
                p <- p + geom_text_repel(data=subset(x,get(paste0(pre.x,"mean", suf[i])) >0),
                                         aes(label=info2),
                                         nudge_y = max.y.axis - unlist(x[get(paste0(pre.x,"mean", suf[i])) >0 , get(paste0(pre.x,"mean", suf[i]))]),
                                         direction="x", vjust=0,
                                         angle = 90,
                                         segment.size  = 0.2,
                                         color='black',
                                         segment.color = "grey50") +
                
                    geom_text_repel(data=subset(x,get(paste0(pre.x,"mean", suf[i])) <0),
                                  aes(label=info2),
                                  nudge_y = min.y.axis - unlist(x[get(paste0(pre.x,"mean", suf[i])) <0 , get(paste0(pre.x,"mean", suf[i]))]),
                                  direction="x", vjust=1,
                                  angle = 90,
                                  segment.size  = 0.2,
                                  color='black',
                                  segment.color = "grey50")

                } else {
                    fp.tp <- lapply(c("yes", "no"), function(t) sapply(c("median", "sd"), function(k) round.f(vector=x[get(null[i-1]) == t & get(null[i]) == "no", get(info.s)], fun=k, digits =2)))

                    labs <- mapply(function(i,j) {
                        paste("Median Info for", i, "positives = " ,paste( j, collapse = ' \u00B1 '))
                    }, c("false", "true"), fp.tp , USE.NAMES=F)
                    
                    p <- p + annotate(geom="text", label= labs, y = c(max.y.axis, max.y.axis - 0.3*max.y.axis), x = min(x$Position)+ 0.01*( max(x$Position) - min(x$Position)), hjust=0, vjust=0)

                }
            }
                                       
            return(p)
    })

    a <- plot_grid(plotlist=p ,ncol=1)
    b= ggdraw(a) + draw_label(gene, y =.97)
    if(!is.null(fisher)){
        b <- b + draw_label(paste0("pmin.fsnp= ", formatC(unique(x[[fisher]]),digits=1 ) ) , x=0.1, y = 0.35)
    }
    return(b)
}


    #' aux function to round mean or sd (or any function output)
    #' @param fun
    #' @param vector to pass function
    #' @param digits to round by
    #' export
    #' return output from fun rounded by digits
    #' round.f()

    round.f <- function(vector, fun, digits){
        tmp <- round(eval(get(fun)(vector)), digits)
        return(tmp)
        }

############ Troubleshooting                                                                     
                                                                     
#' get name of err file for a giving gene
#'
#' Helps with troubleshooting to get the name of a log file given the gene run
#' @param gene
#' @param path path to log files
#' @keywords troubleshooting log files
#' @export
#' @return data table with Gene_id and log file name
#' gene.log()

gene.log <- function(gene,path){
    tmp <- system(paste0("grep ",gene," -l ",path, "/*"), intern=TRUE)
    dt <- data.table(Gene_id=gene, err.file=tmp)
    return(dt)
}

#' get tail of a err file
#'
#' Helps with troubleshooting to get the error resported from  of a log file
#' @param dt data table with Gene_id and full path to err file (column err.file), outtput from gene.log
#' @keywords troubleshooting log files
#' @export
#' @return list with tail of err file
#' tail.log()

tail.log <- function(dt){
    l <- lapply(1:nrow(dt), function(i) system(paste0("tail ",dt$err.file[i]), intern=TRUE))
    names(l) <- dt$Gene_id   
    return(l)
}

#' get reads from bam file mapping a snp
#'
#' Gets reads from bam file mapping a region and filtered by mapq and pheredQ
#' @param bam path to bam file
#' @param snp character vector with chr:pos-pos, chr as formatted in bam file, region to extract reads from
#' @param mapq , whether to select reads meetting a particular mapq threshold, defualts to uniquely mapped as coded in STAR
#' @param Q, whether to select reads when the snp base was called with Q above cut-off, defaults to null
#' @keywords bam reads
#' @export
#' @return data table with selected reads as formatted in bam files
#' bam.r()

bam.r <- function(bam,snp,mapq=255,Q=NULL){

   sam.command <- paste("samtools view ", bam,snp, "| grep ",mapq) 

    tmp <- fread(sam.command, header=F)#, sep=sep, colClasses=c(V4="character", V5="character"))
    ## Don't know how to deal with this now, need to process cigar string and know what to do with the different options, leave this for later.
    ## filter by Q, Q col is always 11 in bam file
    ##QUAL: ASCII of base QUALity plus 33
    ##snp.ps <- as.numeric(gsub("-.*","",gsub(".*:","",snp)))
    ## get snp position for each read
    ##tmp[, snp.pos:= snp.ps-V4]

    return(tmp)
}

#' calculate pair-wise conditional probabilities for fsnps within gene
#'
#' calculates pSNPi== het(hom) | SNPj==het(hom)
#' @param dt data table with GT data with scale 0,1,2 or 0,1,-1,2; NAs allowed. Rows snps, cols CHROM, POS,REF,ALT, samples with sufix "_GT", extra cols allowed
#' @param GT, either het or hom, defaults to het
#' @keywords 
#' @export
#' @return data table ....
#' p.GT()

p.GT <- function(dt, GT="het"){
    samps <- grep("_GT",names(dt),value=T)
    dt[,id:=paste(CHROM,POS,REF,ALT, sep=":")]
    tmp <-as.matrix(dt[,samps,with=F])
    colnames(tmp) <- samps
    rownames(tmp) <- dt$id
    gt <- ifelse(GT=="het",1,0)
    ## create list, each element has the indices to condition by each SNP
    tmp2 <- apply(tmp, 1, function(i) which(abs(i)==gt))
    mat <- matrix(1,nrow(tmp),nrow(tmp), dimnames=list(dt$id, dt$id))
    for(i in 1:nrow(tmp)){
        no.i <- which(!1:nrow(tmp) %in% i)
        mat[no.i,i] <- apply(tmp[no.i,tmp2[[i]], drop=F], 1, function(j) sum(abs(j)==gt, na.rm=T)/sum(!is.na(j)))
    }
    return(mat)
}

#' Get proportion of false positives relative to true positives or consistent positives
#'
#' @param dt data table with stan output, required cols 'col' (below), min.p.fsnp, and cols ending with paste0(".*null",suf)
#' @param col column name from dt to get false positive percentage by, defaults to info
#' @param suf character vector with suffix for cols to compare, first gold standard, second test
#' @param col name of column to get FP(CP), defaults to 'info'
#' @param v vector with cut-off of 'col', when info score, defaults to 0.3,0.4,0.6,0.8,1
#' @param ex whether to exclude info values above ex
#' @param p.fsnp optional number [0,1],cut-off to exclude associations with min.p.fsnp below given value
#' @param pos whether to use CP consistent positive or TP true positive in output table
#' @param xtable, optional to print xtable, defaults to "no", NULL
#' @keywords 
#' @export
#' @return data table with cols: v, FP(CP)% and number of TP(CP), optional print xtable 
#' p.fp()

p.fp <- function(dt, suf=c("",".rna"), col='info', v= c(0.3,seq(0.4,1,0.2)),ex, p.fsnp=10^-8,pos="CP", xtable=NULL) {
    
    dt.i <- dt[info<=ex &  min.p.fsnp>=p.fsnp,]
    tp <- c()
    fp.per <- c()
    gold.st <- grep(paste0("null",suf[1], "$"), names(dt), value=T)
    test <-grep(paste0("null",suf[2],"$"),names(dt), value=T)
        
    for(i in v){
        tmp <- dt.i[get(col)>=i & get(test) =="no" , .N,  by= c(gold.st)]
        temp.tp <- ifelse(length(tmp[get(gold.st) =="no",N]), tmp[get(gold.st) =="no",N], 0)
        tp <- c(tp,temp.tp )     
        fp <- ifelse(length(tmp[get(gold.st) =="yes",N]), tmp[get(gold.st) =="yes",N], 0)    
        fp.per <- c(fp.per, round(100*fp/(temp.tp+fp),0))
    }
    ## create matrix
    FP <- rbind(fp.per,tp)
    neg <- ifelse(pos=="CP","CN","FN")
    ## xtable
    if(!is.null(xtable)) {
        rownames(FP) <- c(paste0(neg,"(%)"), pos)

        addtorow <- list()
        addtorow$pos <- list(0, 0)
        addtorow$command <- c(paste0("& \\multicolumn{5}{c}{",col, "$\\geq$} \\\\\n"),
                              paste("&", paste(v, collapse=" & "),  " \\\\\n"))

        print(xtable(FP), add.to.row = addtorow, include.colnames = FALSE,  booktabs = TRUE, size="small")
    }

    rownames(FP) <- c(paste0(neg,".per"), pos)    
    fp <- data.table(t(FP))
    fp[, eval(col):= v]
    return(fp)


    
}

#' Get proportion of false (consistent) positives relative to true positives or consistent positives combining a range of CIs for RNA and DNA
#'
#' @param dt data table with stan output, required cols 'col' (below), min.p.fsnp, and cols ending with paste0(".*null",suf)
#' @param col name of column to get FP(CP), defaults to 'info'
#' @param v vector with cut-off of 'col', when info score, defaults to 0.3,0.4,0.6,0.8,1
#' @param ex whether to exclude info values above ex
#' @param p.fsnp optional number [0,1],cut-off to exclude associations with min.p.fsnp below given value
#' @param pos whether to use CP consistent positive or TP true positive in output table
#' @param xtable, optional to print xtable, defaults to "no", NULL
#' @param suf.dna suffixes to 'xxxnull' for DNA or gold standard
#' @param suf.rna suffixes to 'xxxnull' for RNA or test
#' @param ci.dna vector with values for gold standard CI, if null suf.dna is used
#' @param ci.rna vector with values for test CI, if null suf.rna is used
#' @param ref name for gold standard, defaults to DNA
#' @param test name for test, defaults to RNA
#' @keywords 
#' @export
#' @return data table with cols: DNA.CI, v, FP(CP)% , number of TP(CP),RNA.CI ,optional print xtable 
#' comb.fp()

comb.fp <- function(dt, col='info', v= c(0.3,seq(0.4,1,0.2)),ex, p.fsnp=10^-8,pos="CP", xtable=NULL, suf.dna, suf.rna, ci.dna=NULL, ci.rna=NULL, ref="DNA", test="RNA"){

    l <- list()
    if(is.null(ci.dna)) ci.dna =suf.dna
    if(is.null(ci.rna)) ci.rna = suf.rna

    for (x in seq_along(suf.rna)) {        
        cp <- lapply(suf.dna, function(i) p.fp(dt, suf=c(i,suf.rna[x]), col=col,v= v,ex=ex, p.fsnp=p.fsnp,pos=pos))

        names(cp) <- as.character(ci.dna)
        cp <- rbindlist(cp, idcol=paste0(ref,".CI"))
        l[[x]] <- cp
    }
    names(l) <- as.character(ci.rna)
    tmp <- rbindlist(l, idcol=paste0(test,".CI"))
    return(tmp)
}

#' Plot proportion of false (consistent) positives relative to true positives or consistent positives combining a range of CIs for RNA and DNA
#'
#' @param cpRNA data table test CI limits (RNA),gold standard CI limits (DNA), %CP, number of true pos, col to assess by, defaults to info, output from com.fp
#' @param col name of column to get FP(CP), defaults to 'info'
#' @param v vector with values of col to add CP/TP in plot, defaults to NULL
#' @param ref col name for gold standard CI limits, defaults to DNA.CI
#' @param test col name for test CI limits, defaults to RNA.CI
#' @param xlab string with name for x axis, defaults to "Info score"
#' @param y yaxis to plot, defaults to CN.per (consistent negative %)
#' @param remove.legend, whether to remove legend, defaults to NULL (with legend)
#' @keywords 
#' @export
#' @return ggplot object
#' plot.fp()

plot.fp <- function(cpRNA,col='info',v=NULL,  ref="DNA.CI", test="RNA.CI", xlab="Info score", y="CN.per", remove.legend=NULL) {

    cpRNA[, eval(test):= paste0(eval(test)," = ", get(test))]

    ## add number of consistent positives for selected values of info for plotting
    if(!is.null(v)){
        cpRNA[, plot.lab:= ifelse(get(col) %in% v, CP,"")]

        ## remove repeated entries to simplify plot
        setkeyv(cpRNA, c(test, col))

        dups <- rle(cpRNA$plot.lab)
        w <- which(dups$length > 1) ## repeated consecutive elements
        l <- dups$length
        for(j in w){
            if(j==1) { ## first element repeated
                cpRNA$plot.lab[2: l[1] ]   <- ""
            } else {
                cpRNA$plot.lab[(sum(l[1:(j-1)]) + 2) : sum(l[1:j]) ]   <- ""
            }                 
        }
    }

    y.lab <- ifelse(y=="CN.per", "Consistent Negatives (%)", "Number of Consistent positives")
    p <- ggplot(cpRNA, aes(get(col), get(y), group=get(ref), color=get(ref), shape=get(ref))) +
        geom_line(linetype=3) +
        geom_point() +
        geom_hline(yintercept=0, linetype=4) +    
        labs(x= paste0(xlab, " \u2265"), y=y.lab)+
        facet_grid(as.formula(paste('.~',test))) + 
        theme_bw()   
    
    if(is.null(remove.legend)){
        p <- p + scale_colour_discrete(name=ref) +
            scale_shape_discrete(name=ref) +
            guides(colour = guide_legend(reverse=T), shape = guide_legend(reverse=T) )
    } else {

        p <- p + theme(
                     legend.text = element_text(color = "white"),
                     legend.title = element_text(color = "white"),
                     legend.key = element_rect(fill = "white")
                 ) + 
            scale_color_discrete(guide = guide_legend(override.aes = list(color = "white", shape = "white")))
    }
      
    if(!is.null(v)){
        p <- p  + geom_text_repel(aes(label=plot.lab))
    }

    return(p)
       ## theme(strip.background = element_rect( fill=NA))
}

#' Get proportion of false (consistent) positives relative to true positives or consistent positives combining a range of CIs for RNA and DNA and plot, wrapper for com.fp and plot.fp
#'
#' @param dt data table with stan output, required cols 'col' (below), min.p.fsnp, and cols ending with paste0(".*null",suf)
#' @param col name of column to get FP(CP), defaults to 'info'
#' @param v vector with cut-off of 'col', when info score, defaults to 0.3,0.4,0.6,0.8,1
#' @param ex whether to exclude info values above ex
#' @param p.fsnp optional number [0,1],cut-off to exclude associations with min.p.fsnp below given value
#' @param pos whether to use CP consistent positive or TP true positive in output table
#' @param suf.dna suffixes to 'xxxnull' for DNA or gold standard
#' @param suf.rna suffixes to 'xxxnull' for RNA or test
#' @param ci.dna vector with values for gold standard CI, if null suf.dna is used
#' @param ci.rna vector with values for test CI, if null suf.rna is used
#' @param ref name for gold standard, defaults to DNA
#' @param test name for test, defaults to RNA
#' @param v2 vector with values of col to add CP/TP in plot, defaults to NULL
#' @param xlab string with name for x axis, defaults to "Info score"
#' @param y yaxis to plot, defaults to CN.per (consistent negative %)
#' @param remove.legend, whether to remove legend, defaults to NULL (with legend)
#' @keywords plot consistent positives RNA DNA
#' @export
#' @return ggplot object
#' comb.plot.fp()

comb.plot.fp <- function(dt, col='info', v= c(0.3,seq(0.4,1,0.2)),ex, p.fsnp=10^-8,pos="CP", suf.dna, suf.rna, ci.dna=NULL, ci.rna=NULL, ref="DNA", test="RNA", v2=NULL, xlab="Info score", y="CN.per", remove.legend=NULL){

    ## get false/consistent positives
    tmp <- comb.fp(dt=dt,col=col,v=v,ex=ex,p.fsnp=p.fsnp,pos=pos,suf.dna=suf.dna,suf.rna=suf.rna,ci.dna=ci.dna, ci.rna=ci.rna,ref=ref,test=test)
    
    ref <- paste0(ref,".CI")
    test <- paste0(test,".CI")

    p <- plot.fp(tmp,col=col,v=v2,  ref=ref, test=test, xlab=xlab, y=y, remove.legend=remove.legend)
    return(p)
}

#' Get snp distance to gene, distance is measured to the closest to start or end
#'
#' @param dt1 data table with Gene_id, snp column, other cols allowed
#' @param snp name for snp column, snp in the form of pos:ref:alt, defaults to "tag"
#' @param dt2 data table with gene_id, start and end 
#' @keywords distance gene snp
#' @export
#' @return dt same as dt1 plus extra column named gene.dist
#' gend.d()

gene.d <- function(dt1, dt2, snp="tag"){
    dt <- merge(dt1, dt2, by.x="Gene_id", by.y="gene_id", all.x=T)
    dt[, pos:=as.numeric(sub(":.*","", get(snp)))]
    dt[, st.dist:= pos-start]
    dt[, end.dist:=pos-end]   
    dt[, gene.dist:=st.dist][abs(st.dist) > abs(end.dist), gene.dist := end.dist]
    dt[pos >= start & pos <= end, gene.dist:=0]
    ## for snsp within gene dist=0
    dt[, c("pos","st.dist", "end.dist", "start", "end", "chrom"):=NULL]
    return(dt)
    
}


#' plots snp gene distance from Btrecase output
#'
#' This function allows you to  compare gene-snp distance
#' @param dt data table with gene.dist column, output from gene.d function
#' @param x1 name of col with  null for condition1
#' @param x2 name of col with  null for condition2
#' @param model character vector with model name to associate with x1 and x2 respectively
#' @param xl character with xlab
#' @keywords Btrecase summary plot  
#' @export
#' @return ggplot object
#' d.plot()

d.plot <- function(dt, x1=NULL, x2=NULL, xl, model){   

    p <- ggplot(dt, aes(gene.dist/1000)) +      
        theme_bw() +
        xlab(xl) +
        ylab(NULL) +
        theme(axis.title = element_text(size=30), axis.text.x = element_text(colour="black", size = 28), legend.title = element_text(size=22), 
              legend.text = element_text(size=18)) +     
        guides(colour=guide_legend(override.aes = list(size = 10)))

    if(!is.null(x1) & !is.null(x2)){
        null.cols <- c(x1, x2)
        
        dt[, Signif:="None"][get(null.cols[1])=="no" & get(null.cols[2])=="no", Signif:="Both"][get(null.cols[1])=="no" & get(null.cols[2])=="yes", Signif:=model[1]][get(null.cols[1])=="yes" & get(null.cols[2])=="no", Signif:=model[2]]
        dt[,Signif:=factor(Signif, levels=c("None", model, "Both"))]
        p <- p + geom_density(aes(colour=Signif, fill=Signif), alpha=0.05)
    } else {

        p <- p + geom_density()
    }
    return(p)
    
}

    
    


    
        

 
