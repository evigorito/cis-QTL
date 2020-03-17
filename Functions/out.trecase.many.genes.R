library(data.table)
library(parallel)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggrepel)
library(gtable)
library(mixtools)
library(grid)
library(hrbrthemes)
library(lemon)

## Functions to help processing stan output

source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

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
    tmp <- rbindlist(lapply(tmp,fread), fill=TRUE)
    return(tmp)
}




#' combine summaries from stan models allowing different tags in 1 model
#'
#' This function allows you to combine summaries from diffrent stan models allowing for different tags in 1 model
#' @param x data table with summary for model 1, same tags as model 2
#' @param y data table with summary for model 2, same tags as model 1
#' @param z data table with summary for model 3, some share some different tags
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
    n <- grep(paste0("se_mean", s[1]), names(all), value=T)
    if(length(grep("log2\\(", n))) {
        no.match.z <- all[is.na(`log2(aFC)_mean.ngt`) , .(Gene_id,tag)] ## tags run in gt/lm
        no.match.y <- all[is.na(`log2(aFC)_mean.gt`) , .(Gene_id,tag)] ## tags run in ngt
    } else {
        m3 <- grep(paste0("log2_aFC_mean.*",s[3]), names(all), value=T)
        m2 <- grep(paste0("log2_aFC_mean.*",s[2]), names(all), value=T)
        no.match.z <- unique(all[is.na(get(m3[1])) , .(Gene_id,tag)]) ## unique tags run in gt/lm, more general form
        no.match.y <- unique(all[is.na(get(m2[1])) , .(Gene_id,tag)]) ## tags run in ngt

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

    x.y.z[,paste0("tag",s[2]):=tag][,paste0("tag", s[3]):=tag][,tag:=NULL]

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



            
    

#' given a tag snp in one set it looks at if at least one of the tagged snps is part of a group in a different set and if so gives the gene, tag in set 1 and tag in set2. 
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
#' This function allows you to compare account for definition of REF allele between a tag SNP and the SNP tagged in the direction of effects
#' @param DT data table with 2 tag cols plus stan output
#' @param le.file full path to legend file for the corresponding chr to look-up for eaf
#' @param s vector with suffixes for tag column names
#' @keywords stan output tag
#' @export
#' @return data table with extra cols with eaf for tags and corrected direction of effect size
#' eaf.aj()

eaf.aj <- function(DT,le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz',s){
    eaf.dif.tags <- apply(DT[same.tag=="no", paste0("tag",s), with=FALSE], 2, snp.eaf, file1=le.file)

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
#' @param path.t path to files with tag files for nogt model, defaults to NULL if same as path.s
#' @param pattern.t pattern for tag files
#' @param lm.sum data table with summary of lm
#' @param gt.sum data table with summary of gt
#' @param s vector with suffixes for tag column names for lm, gt and ngt
#' @param tags.m2 data table with tags run in model 2 (gt)
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
        tmp[, op.dir:='no'][(get(paste0("tag.EAF",s[2]))<0.5 & get(paste0("tag.EAF",s[3])) >0.5) | (get(paste0("tag.EAF",s[2])) >0.5 & get(paste0("tag.EAF",s[3])) <0.5), op.dir:="yes"]
        col.ngt <- grep(paste0("log2.*[^se]_mean", s[3]), names(tmp), value=T)
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
#' @return data table input with new column
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
#' @param s number of row to sample when data table is too long, defaults to 5000
#' @param rx vector with range for x axis
#' @param ry vector with range for y axis
#' @param xl character with xlab
#' @param yl character with ylab
#' @param col when usimg color for symbol based on sig or not sig, provide identifier for first and second column being significant
#' @param title character with title for plot, optional
#' @param title.size size for title
#' @param axis.title size for axis.title 
#' @param axis.text size for axis.text.x and axis.text.y 
#' @param legend.title size for legend title
#' @param legend.text size for legend.text
#' @param legend.symbol size for legend symbols
#' @param point size for geom_points
#' @param colors colors to use for symbols in plot, defaults to c("#999999", "#F0E442", "#0072B2", "#D55E00") for each level of Significance
#' @param shapevar name of variable to use for shape, defaults to NULL
#' @param shape character vector with shape to use, names factor in shapevar, values shape, defaults to NULL
#' @param shape.leg whether to show shape legend, defaults to TRUE
#' @keywords Btrecase summary plot  
#' @export
#' @return ggplot object
#' btrecase.plot()

btrecase.plot <- function(dt, x1, x2, s=5000, rx=NULL,ry=NULL,xl, yl, col=NULL, title=NULL, title.size=12,axis.title=12, axis.text=10, legend.title=12, legend.text=10, legend.symbol=4, point.size=3, colors=c("#999999", "#F0E442", "#0072B2", "#D55E00"), shapevar=NULL, shape=NULL, shape.leg="legend"){
    
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
    if(nrow(x)> s) { ## sample s rows
        x <- x[sample(nrow(x), s),]
    }
    
    new.n <-  sub.grep(names(x), "%", "")
    if(length(new.n)){
        setnames(x, old.n(names(x), "%"), new.n)
        x1 <- sapply(x1, gsub, pattern="%",replacement= "", USE.NAMES=F)
        x2 <- sapply(x2, gsub, pattern="%",replacement= "", USE.NAMES=F)

    }
    if(any(names(x)=="op.dir")){
        ## need to make CI quantiles compatible with estimate
        if(length(new.n)){
            q.n <- new.n[new.n %in% x2]
            x[, paste0(q.n,".comp") := lapply(q.n, function(i) x[[i]])]
            x[op.dir=='yes', paste0(q.n,".comp") := lapply(rev(q.n), function(i) -unlist(x[op.dir=='yes',i,with=F]))]
            x[, (q.n) := lapply(paste0(q.n,".comp"), function(i) x[[i]])]
        }
    }



p <- ggplot(x, aes(get(x1[1]), get(x2[1]))) +      
        geom_errorbarh(aes_string(xmin=x1[2], xmax=x1[3]), linetype=2, colour='grey' ,size=0.1)  +
        geom_errorbar(aes(ymin=get(x2[2]), ymax=get(x2[3])),linetype=2, colour='grey', width=0 ,size=0.1)  +
        theme_bw() +
        xlab(xl) + ylab(yl) +
        theme(axis.title = element_text(size=axis.title), axis.text.x = element_text(colour="black", size = axis.text),axis.text.y = element_text(colour="black", size = axis.text), legend.title = element_text(size=legend.title), 
              legend.text = element_text(size=legend.text)) +
        geom_rug(col=rgb(.7,0,0,alpha=.2)) +
        geom_vline(xintercept=0,color = "blue") +
        geom_hline(yintercept=0,color = "blue") +
        geom_abline(intercept=0,slope=1, linetype=2, size=0.5) 
        ##guides(colour=guide_legend(override.aes = list(size = legend.symbol)))
   
    if(!is.null(col)){
        null.cols <- c(x1[4], x2[4])
        
        x[, Signif:="None"][get(null.cols[1])=="no" & get(null.cols[2])=="no", Signif:="Both"][get(null.cols[1])=="no" & get(null.cols[2])=="yes", Signif:=col[1]][get(null.cols[1])=="yes" & get(null.cols[2])=="no", Signif:=col[2]]
        x[,Signif:=factor(Signif, levels=c("None", col, "Both"))]
    
        names(colors) <- levels(x[["Signif"]])
        man.col <- colors[names(colors) %in% unique(x$Signif)]

        if(is.null(shapevar)){
            p <- p + geom_point(aes(colour=Signif), shape=1, size=point.size) +
                scale_colour_manual(values=man.col) ##+
        } else {
            man.shape <- shape[names(shape) %in% x[[shapevar]]]
            p <- p  + geom_point(aes(colour=Signif, shape=get(shapevar)),
                                 size=point.size) +
                scale_colour_manual(values=man.col) +
                scale_shape_manual(values=man.shape) +
                labs(shape=shapevar)
        }
            ##geom_text_repel(aes(label=ifelse(get(x1[1]) %in% x[Gene_id %in% fpg, get(x1[1])], Gene_id, "")))
    } else {
        if(is.null(shapevar)){
            p <- p + geom_point(shape=1, size=point.size)
        } else {
             man.shape <- shape[names(shape) %in% x[[shapevar]]]
             p <- p  + geom_point(aes(shape=get(shapevar)),
                                 size=point.size) +
                        scale_shape_manual(values=man.shape) +
                        labs(shape=shapevar)
        }
    }
    
#############################       
    ## p <- ggplot(x, aes(get(x1[1]), get(x2[1]), label=Gene_id)) +      
    ##     geom_errorbarh(aes_string(xmin=x1[2], xmax=x1[3]), linetype=2, colour='grey' ,size=0.1)  +
    ##     geom_errorbar(aes(ymin=get(x2[2]), ymax=get(x2[3])),linetype=2, colour='grey', width=0 ,size=0.1)  +
    ##     theme_bw() +
    ##     xlab(xl) + ylab(yl) +
    ##     theme(axis.title = element_text(size=axis.title), axis.text.x = element_text(colour="black", size = axis.text),axis.text.y = element_text(colour="black", size = axis.text), legend.title = element_text(size=legend.title), 
    ##           legend.text = element_text(size=legend.text)) +
    ##     geom_rug(col=rgb(.7,0,0,alpha=.2)) +
    ##     geom_vline(xintercept=0,color = "blue") +
    ##     geom_hline(yintercept=0,color = "blue") +
    ##     geom_abline(intercept=0,slope=1, linetype=2, size=0.5) +
    ##     guides(colour=guide_legend(override.aes = list(size = legend.symbol)))

    ## if(!is.null(col)){
    ##     null.cols <- c(x1[4], x2[4])
        
    ##     x[, Signif:="None"][get(null.cols[1])=="no" & get(null.cols[2])=="no", Signif:="Both"][get(null.cols[1])=="no" & get(null.cols[2])=="yes", Signif:=col[1]][get(null.cols[1])=="yes" & get(null.cols[2])=="no", Signif:=col[2]]
    ##     x[,Signif:=factor(Signif, levels=c("None", col, "Both"))]
    
    ##     names(colors) <- levels(x[["Signif"]])
    ##     man.col <- colors[names(colors) %in% unique(x$Signif)]
        
    ##     p <- p + geom_point(aes(colour=Signif), shape=1, size=point.size) +
    ##         scale_colour_manual(values=man.col) ##+
    ##     ##geom_text_repel(aes(label=ifelse(get(x1[1]) %in% x[Gene_id %in% fpg, get(x1[1])], Gene_id, "")))
    ## } else {   
    ##     p <- p + geom_point(shape=1, size=point.size)
    ## }

    
        
    if(!is.null(title)){
        p <- p + ggtitle(title) +
            theme(plot.title=element_text(size=title.size,hjust = 0.5))
        
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
#' @param ci whether to add ci to plots, NULL is yes
#' @param ci.x  numeric with probability of ci containing true value for gt, defaults to 0.95 (2.5 and 97.5)
#' @param null.x column with null to use if more than one, defualts to NULL
#' @param y data table with stan output for ngt and rna, output of comp.ngt or var.e merged with rna
#' @param ci.y  numeric with probability of ci containing true value for rna, defaults to 0.95 (2.5 and 97.5)
#' @param null.y character with the prefix of the column with null to use if more than one, defualts to NULL
#' @param ci.z  numeric with probability of ci containing true value for rna, defaults to 0.95 (2.5 and 97.5)
#' @param yvar name of variable to plot in y axis, defaults to log2_aFC_mean
#' @param yaxis name for y axis
#' @param xcoord y coordinate for adding text label to plot, defaults to 0.1
#' @param ycoord, y coordinate for adding text label to plot, defaults to 0.2
#' @param colvar name of variable to use for color, defaults to Signif
#' @param colors.x character vector with colors to use, names factor in colvar for x, values colors
#' @param colors.y character vector with colors to use, names factor in colvar for y, values colors
#' @param sizevar name of variable to use for size, defaults to NULL
#' @param size character vector with size to use, names factor in sizevar, values size, defaults to NULL
#' @param size.leg whether to show size legend, defaults to TRUE
#' @param shapevar name of variable to use for shape, defaults to NULL
#' @param shape character vector with shape to use, names factor in shapevar, values shape, defaults to NULL
#' @param shape.leg whether to show shape legend, defaults to TRUE
#' @param z data table with stan output for all combined to indicate comparable snps
#' @param info.s optional to add info score on second plot, give name of column with info
#' @param fisher optional to add lowest pvalue for fisher test of p(het) between RNA and DNA, give name of column with fisher test
#' @param gene character vector with gene ID
#' @param gene.track ggplot object with gene track, defaults to NULL
#' @param rsid file with variant information formatted as ensembl ftp.ensembl.org/pub/grch37/current/variation/gvf/homo_sapiens/homo_sapiens-chr22.gvf.gz for appropiate built and chromosome, defaults to none
#' @param hline, intercept to add hline to indicate significance, defaults to NULL
#' @keywords stan plot gene-snps compare
#' @export
#' @return ggplot object with gt and rna plots only
#' gene.plot2b()

gene.plot2b <- function(x,ci=NULL, ci.x=NULL, null.x=NULL, y, ci.y=NULL, null.y=NULL, ci.z=NULL, yvar="log2_aFC_mean", yaxis="eQTL-effect",xcoord=0.1, ycoord=0.2, colvar="Signif", colors.x= setNames(c("goldenrod4", "steelblue4","lightgoldenrod1", "steelblue1"),  c("No, Both", "Yes, Both", "No, obs-GT", "Yes, obs-GT")) ,
                        colors.y= setNames(c("goldenrod4", "steelblue4","lightgoldenrod1", "steelblue1"), c("No, Both",  "Yes, Both" , paste0(c('No', 'Yes'), ", hidden-GT"))) ,
                        sizevar=NULL, size=NULL, size.leg="legend",
                        shapevar=NULL, shape=NULL, shape.leg="legend",
                        z,gene,info.s=NULL,
                        fisher=NULL,
                        gene.track=NULL,
                        rsid=NULL,
                        hline=NULL){
    ## select entries for gene
    x <- x[Gene_id==gene,]
    y <- y[Gene_id==gene,]
    ## create col for position for gt and ngt
    x[, Position:=as.numeric(gsub(":.*","", tag))/10^6]
    ## get col names for "tag" in y DT
    tag.y <- grep("tag",names(y),value=T)
    ## use first one
    y[, Position:=as.numeric(gsub(":.*","", get(tag.y[1])))/10^6]
    ## order `log2(aFC)_null.gt` to make consistent coloured plots when only "yes" is present
    if(is.null(null.x)){
        null <- grep("null", names(x),value=T)
    } else {
        null <- null.x
    }
       
    ## merge x with z to get common snps, null column for x from z
    x <- merge(x, z[,c('Gene_id', 'tag.gt', 'tag.ngt', paste0(null,'.gt')), with=F], by.x=c('Gene_id','tag'), by.y=c('Gene_id', 'tag.gt') , all.x=T)
    ## for snps in z (comparable SNPs) change position to nogt tag to ease comparison in plots, but only if the tag.ngt is in y (some were removed by cis-window threshold)

    tags.ngt.out.cis <- x[!is.na(tag.ngt) & !tag.ngt %in% y$tag,tag.ngt]
    x[!is.na(tag.ngt) & tag.ngt %in% y$tag, Position:=as.numeric(gsub(":.*","", tag.ngt))/10^6]
    
    ## need to make NAs tags.ngt.out.cis as they wont have a counterpart in y
    x[tag.ngt %in% tags.ngt.out.cis, paste0(null, ".gt") := NA]
    x[is.na(get(paste0(null, '.gt'))) , paste0(null, ".gt") := paste0(get(null), '.gt')]
    x[, paste0(null, ".gt"):=factor(get(paste0(null, ".gt")), levels=c("yes", "no", "yes.gt", "no.gt"))]
    x[,(null):= get(paste0(null, ".gt"))]

    x[, paste0(null, ".gt"):=NULL]

   

    ## get suffixes in y
    if(is.null(null.y)){
        suf <- gsub('.*_null', "", grep('.*_null', names(y), value=T))
        ## get prefixes for null in y
        pre <- unique(mapply(function(x,y) {gsub(y,"",x)}, grep('.*_null', names(y), value=T),suf,USE.NAMES=F))
    } else{
        suf <- gsub(null.y, "", grep(null.y, names(y), value=T))
        pre <- null.y
    }

    ## merge y with z, keep eQTL effect corrected from z to make it compatible  with GT in plot
    y <- merge(y,z[,c('Gene_id', 'tag.ngt', "tag.gt", "op.dir",  unlist(lapply(c(pre,"log2_aFC_mean"), function(i)  paste0(i ,c('.ngt','.rna'))))), with=F], by.x=c('Gene_id','tag'), by.y=c('Gene_id', 'tag.ngt') , all.x=T, suffixes=c('',2))

   

    ## replace NAs from z with value from relevant y col
    for(i in suf){
        tmp <- paste0(pre,i,2)
        y[is.na(get(tmp)), (tmp) := paste0(get(paste0(pre, i)),i)]
        lev <- paste0(c("yes", "no"), rep(c("", i), each=2))
        names(lev) <- paste0(c("No", "Yes"), rep(c(", Both", ", hidden-GT"), each=2))
        y[ , (tmp):=factor(get(tmp), levels=paste0(c("yes", "no"), rep(c("", i), each=2)))]
        ## Add Signif column for simplicity
        for( n in seq_along(lev)){          
            y[get(tmp) == lev[n] , paste0("Signif", i):= names(lev[n])]
        }
    }

        
    ## make plots with same x axis and same y axis
    min.pos <- min(c(x$Position, y$Position))
    max.pos <- max(c(x$Position, y$Position))


    ## get prefixes for x and y
   
    pre.x <- gsub('_se_mean', "", grep('.*_se_mean', names(x), value=T))
    pre.y <- unique(gsub('_se_mean.*', "",grep('.*_se_mean', names(y), value=T)))

    ## get ci limits
    p.x <- ifelse(is.null(ci.x), 0.95, ci.x)
    p.y <- ifelse(is.null(ci.y), 0.95, ci.y)
    p.z <- ifelse(is.null(ci.z), 0.95, ci.z)
    p.yz <- c(p.y,p.z)

    if(is.null(ci)){
        min.y.axis <- min(c(x[[paste0(pre.x,'_',(1-p.x)*50,'%')]], sapply(1:length(suf), function(i) y[[paste0(pre.y,  '_',(1-p.yz[i])*50,'%' , suf[i])]])))
        max.y.axis <- max(c(x[[paste0(pre.x,'_',(p.x+(1-p.x)/2)*100, '%')]], sapply(1:length(suf), function(i) y[[paste0(pre.y, '_',(p.yz[i] +(1-p.yz[i])/2)*100,'%' , suf[i])]])))
    } else {
        max.y.axis <- max(c(x[[yvar]], y[[paste0(yvar, suf[1])]], y[[paste0(yvar, suf[2])]])) + 0.2*max(c(x[[yvar]], y[[paste0(yvar, suf[1])]], y[[paste0(yvar, suf[2])]]))
        min.y.axis <- 0
    }
    
    ## make Signif columns for simplicity
    x[, Signif:="No, Both"]
    x[get(pre) == "yes.gt", Signif:="No, obs-GT"][get(pre)=="no.gt", Signif:="Yes, obs-GT"][get(pre)=="no",Signif:="Yes, Both"]
    ##colors <- c("goldenrod4", "steelblue4","lightgoldenrod1", "steelblue1")
    ## colors <- c("red", "blue4", "lightpink3", "skyblue4")
    ## names(colors) <- c("No, Both", "Yes, Both", "No, obs-GT", "Yes, obs-GT")
    man.col <- colors.x[names(colors.x) %in% x[[colvar]]]

    gene.name <- unique(z[Gene_id==gene,Gene_name])

    ## add to plot "GT"
    grob <- grobTree(textGrob("Observed-GT",x=xcoord,  y=ycoord, hjust=0),
                     gp=gpar(fontsize=10))

    ## Add rsid to top cis-SNP, otherwise make col with NA as wont interfere
    if(!is.null(rsid)){
        tmpx <- x[get(null) == "no",][, abs.mean:=abs(log2_aFC_mean)]
        if(nrow(tmpx)){
            if(yvar == "log2_aFC_mean"){
                ## need to choose rsid as the sig snps with max abs(yvar) in this case
                setorder(tmpx, -abs.mean)
            } else {
                ## need to order first by yvar but if more than one with top value then I choose abs.mean
                setorderv(tmpx, c(yvar, "abs.mean"), order=rep(-1,2))
            }
            top.tag <- tmpx[1, tag]
            rs <- get.rsid(rsid, snps=top.tag)
            x <- merge(x, rs, by.x="tag", by.y="SNP", all.x=T)
            ## if top tag is repeated (due to comparison with no gt, make sure only one is selected
            rowpos <- x[!is.na(rsid) , which=T]
            if(length(rowpos)>1) x[rowpos[-1], rsid:=""]
            ## select the rsSNP once for labelling
            rows <- x[is.na(rsid) , which=T]
            x[ rows, rsid:=""]
        } else {
            x[, rsid:=""]       
        } 
           
    } else {
        x[, rsid:=""]
    }
    
    ## x[is.na(rsid), rsid:=""]
    setkey(x,rsid)
    lev.rsnp <- levels(as.factor(x$rsid))

    if(nlevels(as.factor(x$rsid))== 2){
        shape <- rep(1,2)
        shape[which(lev.rsnp != "")] <- 23
        size.p <- rep(1.5,2)
        size.p[which(lev.rsnp != "")] <- 3
        names(shape) <- names(size.p) <- lev.rsnp
    } else {
        shape <- 1
        size.p <- 1.5
        names(shape)  <- names(size.p) <- lev.rsnp
    }

         
    leg.size=9
    axis.txt=8
    symb.size=2

    
    p <- ggplot(x, aes(x=Position, get(yvar) , color=Signif, shape=rsid,
                       size=rsid))  +
        geom_hline(yintercept=0)
    if(is.null(ci)){
        p <- p+ geom_errorbar(data=x, aes(ymin=get(paste0(pre.x,'_',(1-p.x)*50,'%')), ymax=get(paste0(pre.x,'_', (p.x+(1-p.x)/2)*100, '%'))), color="grey90",linetype="dashed" , width=0.001, size=0.2)
    }
    
    p <- p + xlim(min.pos-1/10^5,max.pos+1/10^5) +
                                        #labs(color="Significant")+
        geom_point() +
        geom_text_repel(data=x, aes(label=rsid)) +
        scale_color_manual(name="Signif & Tested", values=man.col) +
        scale_shape_manual(values=shape) +
        scale_size_manual(values=size.p) +
        annotation_custom(grob) +
        theme(axis.title=element_blank(), 
              axis.text.x=element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y=element_text(size=axis.txt),
              legend.title=element_text(size=leg.size),
              legend.text=element_text(size=leg.size)) +
        guides(colour=guide_legend(override.aes= list(size=symb.size, shape=1)),
               shape=FALSE,
               size=F)

    cond <- c("DNA pha-SNPs", "Hidden-GT")
    ##names(colors) <- c("No, Both",  "Yes, Both" , paste0(c('No', 'Yes'), ", hidden-GT")) 
    
    col.cond <- lapply(suf, function(i) colors.y[names(colors.y) %in% y[[paste0("Signif", i)]]])

    ## add cond to plot 
    grob <- lapply(cond, function(i) grobTree(textGrob(i, x=xcoord,  y=ycoord, hjust=0),
                                              gp=gpar(fontsize=10)))


    if(!is.null(rsid)) {
        top.tag <- lapply(suf, function(i){
            g1 <- paste0(null,i,2)
            g2 <- paste0( "log2_aFC_mean", i,2)
            tmpy <- y[get(g1)=="no"  ][,"abs.mean":=abs(get(g2))]
            if(yvar == "log2_aFC_mean"){
                ## need to choose rsid as the sig snps with max abs(yvar) in this case
                setorder(tmpy, -abs.mean)
            } else {
                setorderv(tmpy, c(paste0(yvar,i), "abs.mean"), order=rep(-1,2))
            }            
            return(tmpy[1, tag])
            ##return(unique(y[get(g1) == "no", .(abs(get(g2)), tag)][V1 == max(V1),tag][1]))
        })
        if(length(unlist(top.tag))){
            rs <- get.rsid(rsid,snps= unique(unlist(top.tag)))
        }
       
        tmp <- mapply(function(i,j) {
            y[, paste0("rsid", j):=""]
            if(length(i)){
                y[tag==i,  paste0("rsid", j):=rs[SNP==i,rsid]]
            }},          
            i=top.tag,
            j=suf, SIMPLIFY=F)

    } else {

        tmp <- lapply(suf, function(i) y[, paste0("rsid", i):=""])
    }


    lev.snp <- lapply(suf, function(i) levels(as.factor(y[[paste0("rsid",i)]])))

    shape <- lapply(lev.snp, function(i) {
        if(length(i) ==2){
            shape <- rep(1,2)
            shape[which(i != "")] <- 23
        } else {
            shape <- 1
        }
        names(shape) <- i
        return(shape)
    })

    size.p <- lapply(lev.snp, function(i) {
        if(length(i) ==2){
            shape <- rep(1.5,2)
            shape[which(i != "")] <- 3
        } else {
            shape <- 1.5
        }
        names(shape) <- i
        return(shape)
    })
                                                          

    p23 <- lapply(1:length(suf),  function(i) {
 
        ## make eQTL and CIs signs  compatible with GT tags
        y[get(paste0("Signif", suf[i])) == "Yes", paste0(pre.y,"_mean", suf[i]):= get(paste0(pre.y,"_mean", suf[i], "2"))]
        y[get(paste0("Signif", suf[i])) == "Yes" & sign(get(paste0(pre.y, '_',(1-p.yz[i])*50,'%', suf[i]))) != sign(get(paste0(pre.y,"_mean", suf[i]))), paste0(pre.y, '_',(1-p.yz[i])*50,'%', suf[i]) := - get(paste0(pre.y, '_',(1-p.yz[i])*50,'%', suf[i]))]

        y[get(paste0("Signif", suf[i])) == "Yes" & sign(get(paste0(pre.y, '_',(1-p.yz[i])*50,'%', suf[i]))) != sign(get(paste0(pre.y,"_mean", suf[i]))), paste0(pre.y,'_',(p.yz[i] +(1-p.yz[i])/2)*100,'%', suf[i]) := -get(paste0(pre.y,'_',(p.yz[i] +(1-p.yz[i])/2)*100,'%', suf[i]))]

        p <- ggplot(y, aes(x=Position, get(paste0(yvar, suf[i])),
                           color=get(paste0("Signif", suf[i])),
                           shape=get(paste0("rsid", suf[i])),
                           size=get(paste0("rsid", suf[i])))) 
        if(is.null(ci)){
            
            p <- p + geom_errorbar(data=y, aes(ymin=get(paste0(pre.y, '_',(1-p.yz[i])*50,'%', suf[i])), ymax=get(paste0(pre.y,'_',(p.yz[i] +(1-p.yz[i])/2)*100,'%', suf[i]))), color="grey90",linetype="dashed" ,width=0.001, size=0.2)
        }
        
        p <- p +  xlim(min.pos-10/10^6,max.pos+10/10^6) +
            ylim(min.y.axis-0.05, max.y.axis+0.05) +
            #labs(color="Significant") +
            ##labs(x="Position (Mb)", y=paste0("log2(aFC): ",cond[i]), color="Significant") +
            geom_point() +
            scale_color_manual(name="Signif & Tested", values=col.cond[[i]]) +
            scale_shape_manual(values=shape[[i]]) +
            scale_size_manual(values=size.p[[i]]) +
            
            geom_text_repel(data=y, aes(label=get(paste0("rsid", suf[i])))) +
            ##ggtitle("")
            theme(axis.title=element_blank(), 
                  axis.ticks.x = element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_text(size=axis.txt),
                  legend.title=element_blank(),
                  legend.text=element_text(size=leg.size)) +
            guides(colour=guide_legend(override.aes= list(size=symb.size, shape=1)),
                   shape=F, size=F)
                  
        if(i==1 & !is.null(info.s)){
                    labs <- paste("Median info  =" , round(mean(y[['info.ngt']]),1),'\u00B1', round(sd(y[['info.ngt']]),1))
                    p <- p + annotate(geom="text", label= labs, y = c(max.y.axis), x = min(x$Position)+ 0.01*( max(x$Position) - min(x$Position)), hjust=0, vjust=0)

        }
        if(i==2 & !is.null(fisher)){
                    labs <- paste("min p-value phasing-SNP  =" , formatC(unique(y[[fisher]]),digits=1 ))
                    p <- p + annotate(geom="text", label= labs, y = c(max.y.axis), x = min(x$Position)+ 0.01*( max(x$Position) - min(x$Position)), hjust=0, vjust=0)

        }

        p <- p + annotation_custom(grob[[i]])
        return(p)
    })

    if(is.null(hline)) hline=0
    p <- p + geom_hline(yintercept=hline, linetype="dashed")
    p23 <- lapply(p23, function(i) i+ geom_hline(yintercept=hline, linetype="dashed"))
    
    

    if(is.null(gene.track)){
        a <- plot_grid(plotlist= c(list(p,p23[[2]])), ncol=1)
    } else {
        p1 <- ggplot_gtable(ggplot_build(p))
        p2 <- ggplot_gtable(ggplot_build(p23[[2]]))
        p3 <- ggplot_gtable(ggplot_build(gene.track))
        
        p3$widths <- p1$widths
        p2$widths <- p1$widths

        title.grob <- textGrob(gene.name, 
                   gp=gpar(fontsize=10))

        y.grob <- textGrob(yaxis, 
                   gp=gpar( fontsize=11), y = 0.65, rot=90)

        ##x.grob <- textGrob("Position",
                          ## gp=gpar( fontsize=11), y=1)
        
        ##plot1.2 <- plot_grid(p1,p2, ncol=1 )
        ##top <- arrangeGrob(plot1.2, )

        
        plots <- plot_grid(p1,p2, p3, ncol=1,
                       rel_heights=c(1,1,1))

        a <- plot_grid(arrangeGrob(plots,
                         left = y.grob,
                         ##bottom = x.grob,
                         top=title.grob))
        
       
    } 
    
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
#' @param mod character vector with model name to associate with x1 and x2 respectively
#' @param xl character with xlab
#' @keywords Btrecase summary plot  
#' @export
#' @return ggplot object
#' d.plot()

d.plot <- function(dt, x1=NULL, x2=NULL, xl, mod){   

    p <- ggplot(dt, aes(gene.dist/1000)) +      
        theme_bw() +
        xlab(xl) +
        ylab(NULL) +
        theme(axis.title = element_text(size=30), axis.text.x = element_text(colour="black", size = 28), legend.title = element_text(size=22), 
              legend.text = element_text(size=18)) +     
        guides(colour=guide_legend(override.aes = list(size = 10)))

    if(!is.null(x1) & !is.null(x2)){
        null.cols <- c(x1, x2)
        
        dt[, Signif:="None"][get(null.cols[1])=="no" & get(null.cols[2])=="no", Signif:="Both"][get(null.cols[1])=="no" & get(null.cols[2])=="yes", Signif:=mod[1]][get(null.cols[1])=="yes" & get(null.cols[2])=="no", Signif:=mod[2]]
        dt[,Signif:=factor(Signif, levels=c("None", mod, "Both"))]
        p <- p + geom_density(aes(colour=Signif, fill=Signif), alpha=0.05)
    } else {

        p <- p + geom_density()
    }
    return(p)
    
}

#' add CI.ov column to btrecase output testing whether 2CIs overlap
#'
#' Codes the column as "yes" or "no" for overlaps
#' @param dt data table with btrecase output
#' @param ci limist for ci, defaults to .95 if null
#' @param s suffixes for CI columns
#' @keywords Btrecase CI overlap
#' @export
#' @return data table input with new column
#' add.ov()

add.ov <- function(dt, ci=NULL, s){

    ## get names of CI columns by suffix
    p.x <- ifelse(is.null(ci), 0.95, ci)
    ci.lim <- c((1-p.x)*50, (p.x+(1-p.x)/2)*100)
    ci.col <- lapply(s, function(i) sapply(ci.lim, function(j) grep(paste0(j, "%", i, "$"), names(dt), value=T)))
    names(ci.col) <- s

    ## compare
    dt[, CI.ov:="yes"][( get(ci.col[[1]][1]) < get(ci.col[[2]][1]) & get(ci.col[[1]][2]) < get(ci.col[[2]][1]) ) , CI.ov:="no"]
    dt[( get(ci.col[[1]][1]) > get(ci.col[[2]][2]) & get(ci.col[[1]][2]) > get(ci.col[[2]][2]) ) , CI.ov:="no"]

    return(dt)

}    
    
#' calculate FP percentage or biassed estimate percetange from  btrecase output 
#'
#' Codes FP and BP as num% 
#' @param dt data table with btrecase output
#' @param col name of col to calculate FP or BE from (Signif, CI.ov based, assumes Signif has "Both" and CI.ov has "no" and "yes" )
#' @param s suffix for numerator and denominator to calculate FP, defaults to "no" for CI.ov no overlap for BE
#' @param out either FP or BE
#' @param r value to round FP and or BE
#' @keywords Btrecase FP BE
#' @export
#' @return character for FP% or BE%
#' add.ov()

fp.be <- function(dt, col, s="no", out=c("FP", "BE"), r=0){
    
    if(out=="FP"){
        num <- nrow(dt[ get(col) == s[1],])*100
        den <- nrow(dt[ get(col) == s[2] | get(col) =="Both",])
    }
    if(out == "BE"){
        num <- nrow(dt[ get(col) == s,])*100
        den <- nrow(dt)
    }
    return(paste0(round(num/den, r), "%"))
}



#' creates geom for adding table into a facet plot in ggplot
#'
#' https://stackoverflow.com/questions/49475201/adding-tables-to-ggplot2-with-facet-wrap-in-r?noredirect=1&lq=1
#' https://stackoverflow.com/questions/34625165/adding-table-to-ggplot-with-facets (didnt work)
#' @param 
#' @keywords facet table
#' @export
#' @return geom
#' geom_custom()


#######################################################################################
GeomCustom <- ggproto(
  "GeomCustom",
  Geom,
  setup_data = function(self, data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    data
  },

  draw_group = function(data, panel_scales, coord) {
      ## print(data$x)
      ## print(data$y)
      vp <- grid::viewport(x=unique(data$x), y=unique(data$y))
      ##print(data$grob[[1]])
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_custom", g)
  },

  required_aes = c("grob","x","y")

)

geom_custom <-  function(mapping = NULL,
           data = NULL,
           stat = "identity",
           position = "identity",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = FALSE,
           ...) {
    layer(
      geom = GeomCustom,
      mapping = mapping,
      data = data,
      stat = stat,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(na.rm = na.rm, ...)
    )
}

#' calculates p=n/m for each SNP from tot.ase.counts output
#'
#' get n/m per snp per sample
#' @param mat matrix with rownames sample name and cols snp.n and snp.m for counts for allele in hap2 from phaser (n) and total ASE counts (m), output from tot.ase.counts
#' @keywords Btrecase AI p
#' @export
#' @return matrix with colnames samples and rownames snps; values p
#' p.snp()

p.snp <- function(mat){
    ## get m, total ase counts
    m.ase <- nm.snp(mat)
    

    ## get p = n/m
    n.ase <-nm.snp(mat, "n")
    
    p.ase <- n.ase/m.ase
    colnames(p.ase) <- paste0(colnames(p.ase), "_p")
    return(p.ase)
}

#' get  haplotype2  or total ASE reads for each SNP from tot.ase.counts output
#'
#' get n or m reads per snp per sample
#' @param mat matrix with rownames sample name and cols snp.n and snp.m for counts for allele in hap2 as in phaser (n) and total ASE counts (m), output from tot.ase.counts
#' @param reads whether to get hap (n) or total (m) ASE reads
#' @keywords Btrecase AI ASE
#' @export
#' @return matrix with colnames samples and rownames snps; values alt hap or total ase reads
#' nm.snp()

nm.snp <- function(mat, reads=c("m","n")){
    ## get m, total ase counts
    reads <- match.arg(reads)
    ase <- t(mat[, grep(reads, colnames(mat), value=T)])
    rownames(ase) <- sub(paste0("\\.",reads), "", rownames(ase))
    if(reads == "m"){
        colnames(ase) <- paste0(colnames(ase), "_ASE")
    }
    
    return(ase)
}
  
#' plot stan inputs when no GT for QC purposes
#'
#' plots negative binomial and beta binomial inputs for a gene-snp pair
#' @param stan.in stan input prepared by stan.trecase.rna.noGT.eff2 (stan.eff.R)
#' @param case counts for ASE prepared by tot.ase_counts (stan.eff.R)
#' @param tag rSNP run in model to add to plot, defaults to NULL
#' @keywords Btrecase AI ASE
#' @export
#' @return plot stan inputs by genotype
#' plot.in()

plot.in <- function(stan.in, c.ase, tag=NULL){

    ## get samples with ASE info
    #ncase=c.ase[names(stan.in$ase$m),]

    ## get observed AI per sample and estimated refrence panel bias
    obAI <- beta.in(stan.in, c.ase, ASE="observed")
    estAI <- beta.in(stan.in, c.ase, ASE="estimated")

    obAI[, AI:= "observed"]
    estAI[, AI:= "estimated"]

    ## combine observed and estimated

    both <- rbind(obAI, estAI)

    p1 <- ggplot(both, aes(g,propxp, color = AI)) + geom_jitter(shape=1) + geom_hline(yintercept=0.5, linetype="dashed", color = "red") +  ggtitle("BB noGT")+ expand_limits(y=0)

    nb.g <- sapply(stan.in$NB$p.g,function(j) sum(as.numeric(names(j))*j))
    nb.ngt <- data.table(g=nb.g, counts.log=log(unlist(stan.in$NB$counts)))

    ## remove 0 counts

    nb.ngt <- nb.ngt[is.finite(rowSums(nb.ngt)),]

    ## print(length(nb.g))
    
    print(summary(lm(g~counts.log, data=nb.ngt))$coefficients)

    p2 <- ggplot(nb.ngt, aes(x=g,y=counts.log)) + geom_point(shape=1) + geom_smooth(method='lm')+  ggtitle("NB noGT")+ expand_limits(y=0)

    if(!is.null(tag)){
        p2 <- p2 + draw_label(tag, x=0 , y=1.5,hjust=0, size=10)
    }
    

    return(plot_grid(p1,p2))
}

#' helper function to prepare stan BB input when no GT for QC purposes
#'
#' beta binomial inputs for a gene-snp pair
#' @param stan.in stan input prepared by stan.trecase.rna.noGT.eff2 (stan.eff.R)
#' @param case counts for ASE prepared by tot.ase_counts (stan.eff.R)
#' @param ASE. whether to report results for observed ASE per haplotype or estimated allelic imbalance
#' @keywords Btrecase AI ASE
#' @export
#' @return data table with allelic expression for an haplotype for all possible genotypes of rSNP
#' beta.in()

beta.in <- function(stan.in, c.ase, ASE=c("observed", "estimated")){
    ## get samples with ASE info
    ncase=c.ase[names(stan.in$ase$m),]
    ## get AI per sample
    if(ASE == "observed"){
        prop <- lapply(1:length(stan.in$ase$m), function(i) stan.in$ase$n[[i]]/stan.in$ase$m[i])
    } else {
        prop <- stan.in$ase$ai0

    }
    ## correct for gt==-1, prop==1-p
    prop.c <- lapply(1:length(prop), function(i) {
        w1 <- which(stan.in$ase$g[[i]]==-1)
        if(length(w1)){
            prop[[i]][w1] <- 1-prop[[i]][w1]
        }
    
        names(prop[[i]]) <- abs(stan.in$ase$g[[i]])    
        return(prop[[i]])
    })

    ## proportion * p(prop)
    prop.p <- mapply('*', prop.c,stan.in$ase$p, SIMPLIFY=F)
    if(is.numeric(prop.p)){

        prop.het <- prop.p[which(names(prop.p)==1)]
        prop.hom <- prop.p[which(names(prop.p)!=1)]

    } else {
        
        names(prop.p) <- names(stan.in$ase$n)

        prop.het <- sapply(1:length(prop.p), function(i) sum(prop.p[[i]][names(prop.p[[i]])==1])/sum(stan.in$ase$p[[i]][names(prop.p[[i]]) == 1]))
        
        prop.hom <- sapply(1:length(prop.p), function(i) sum(prop.p[[i]][names(prop.p[[i]])==0 | names(prop.p[[i]])==2])/sum(stan.in$ase$p[[i]][names(prop.p[[i]]) == 0 | names(prop.p[[i]])==2]))

    }

    ## remove NAs
    prop.het <- prop.het[!is.na(prop.het)]
    prop.hom <- prop.hom[!is.na(prop.hom)]
    
    prop.nog <- data.table(g=rep(c("hom","het"), times=c(length(prop.hom), length(prop.het))), propxp=c(sort(prop.hom),sort(prop.het)))

    return(prop.nog)
}


    
#' helper function to format output for report for QC purposes
#'
#' give info about inputs for stan model and print plots
#' @param gene input gene name 
#' @param stan.in stan input prepared by stan.trecase.rna.noGT.eff2 (stan.eff.R)
#' @param case counts for ASE prepared by tot.ase_counts (stan.eff.R)
#' @param gtex.clean object with gtex and bayesian estimates
#' @param clean.sum object with bayesian estimates
#' @param skin.type whether to use normal or psoriasis, defaults to both
#' @keywords Btrecase AI ASE summary
#' @export
#' @return plots made by plot.in
#' info.in()

info.in <- function(gene, stan.in, case, gtex.clean, clean.sum, skin.type){

    ## total counts vs haplotypic counts, same for all snps, get first one

    m <- stan.in[[1]]$ase$m
    tot <- unlist(stan.in[[1]]$NB$counts)
    dt <- data.table(Total_counts=tot[names(m)], Mapped_counts=m)

    print(ggplot(dt, aes(x=Total_counts, y=Mapped_counts)) + geom_point(shape=1))
      
    w <- which(names(stan.in) %in% unique(clean.sum[Gene_name==gene & !is.na(fSNP.id) & skin==skin.type ,tag]))

    for( i in w){
        
        
        print(clean.sum[Gene_name==gene & tag == names(stan.in)[i] & skin == skin.type, .(Gene_name,tag,fSNP.id, log2_aFC_mean, tag.fsnp.op, log2_aFC_null)])
       
        print(plot.in(stan.in[[i]], case, names(stan.in)[i]))
    }
    print("Get fSNPs in gtex")
    
    print(gtex.clean[Gene_name==gene & fSNP=="yes" & skin=='normal_skin' ,.(SNP, log2_fSNP, slope.adj, Signif)])


    print("Get proportion of hets per fSNP")
    print(apply(case[,grep(".m", colnames(case)), drop=F], 2 , function(i) sum(i>0)/nrow(case) ))
}



                         

#' Fromat stan output to correct for multiple testing
#'
#' add columns based on a rejection region of choice and % of posterior out of that region of choice to define significant and no significant calls, based on normal approximation to the posterior
#' @param a numeric vector with one or more rejection values
#' @param b stan summary, required cols are log2_aFC_mean and log2_aFC_sd
#' @param cols if names arent log2_aFC_mean and log2_aFC_sd, defaults to NULL
#' @param c % of posterior out of region to consider significant, defaults to 0.9
#' @keywords multiple testing Bayesian
#' @export
#' @return stan output with 2 new columns null.rej="yes"/"no" based on significance and rej.level with the region threshold. When multiple rejection values the output is in long format and post.out with % of posterior out of rejection zone
#' rej.recode()

rej.recode <- function(a, b, cols=NULL, c=0.9){
    l <- lapply(a, function(i) {
        d <- copy(b)
        if(is.null(cols)){
            d[ abs(log2_aFC_mean) > i, post.out:=pnorm(i,abs(log2_aFC_mean), log2_aFC_sd, lower.tail=FALSE)]
            d[ abs(log2_aFC_mean) <= i, post.out:=0 ]
        } else {
            d[ abs(get(cols[1])) > i, post.out:=pnorm(i,abs(get(cols[1])), get(cols[2]), lower.tail=FALSE)]
            d[ abs(get(cols[1])) <= i, post.out:=0 ]
        }
        
        d[,null.rej:="yes"]        
        d[post.out >= c, null.rej:="no"]
        d[, rej.level:=i]
        return(d)
    })
    dt <- rbindlist(l)
    return(dt)
}


#' Compare significance based on mult testing with CI
#'
#' @param a numeric vector with one or more rejection values
#' @param b stan summary, required cols are rej.null and log2_aFC_null, output from rej.recode
#' @keywords multiple testing Bayesian comparison
#' @export
#' @return data table with comparison by rejection threshold
#' ci.mult()

ci.mult <- function(a, b){
    tmp <- rbindlist(lapply(a, function(i) {
        dt <- b[rej.level == i,]
        dt <- dt[,.N,.(null.rej, log2_aFC_null)]
        dt[,rej.level:=i]
    }))
    return(tmp)
}

    
#' Extract reference panel hap estimates from stan input with GT for each rsnp-gene combination
#'
#' @param inp list with inputs made by in.neg.beta.prob.eff2 (coded in stan.eff.R functions)
#' @param genes when names(inp) are not gene names only (prefix,etc) provide a character vector with the names of the genes that correspond to each element of enp (same order), defualts to NULL
#' @keywords reference panel bias QC
#' @export
#' @return data table with fixed effect for ref panel bias estimate, random effect, SD, tag (SNP) id, Gene_id
#' rpb.hap()

rpb.hap <- function(inp, genes=NULL){
    ## Collect ref panel bias main effect and sd from input by gene_id, snp and ind
    a0 <- list()
    sd0 <- list()
    snp <- list()
    ind <- list()

    for(n in names(inp)){ ## gene level
        if(is.null(genes)) {
            g <- n
        } else {
            g <- genes[which(names(inp) == n)]
        }
        ## extract inputs
        l <- lapply(inp[[n]], '[[', "ai0")
        a0[[g]] <- Reduce ("c",l )
        snp[[g]] <- sapply(l, length)
        
        ########### working here ##############
        

    sig <- ase[Gene_id==g & log2_aFC_null.refbias=="no" & model=="ase", SNP]
    
    names(l) <- sig
    snp.sig[[g]] <- sapply(l, length) 
        a0.sig[[g]] <- Reduce("c", l)



        
    l <- lapply(inp[[n]][!names(inp[[n]])%in%sig], '[[', "ai0")
    names(l) <- names(inp[[n]])[!names(inp[[n]])%in%sig]
    snp.nsig[[g]] <- sapply(l, length)
    a0.nosig[[g]] <- Reduce("c", l)
    sd0.sig[[g]] <- Reduce("c", lapply(inp[[n]][names(inp[[n]])%in%sig], '[[', "sdai0"))
    sd0.nosig[[g]] <- Reduce("c", lapply(inp[[n]][!names(inp[[n]])%in%sig], '[[', "sdai0"))
}

}

#' Extract inputs for a gene-snp association from GT to make plots for QC purposes
#'
#' @param inp stan input for a gene-snp prepared by stan.trecase.eff2 (stan.eff.R)
#' @keywords reference panel bias QC
#' @export
#' @return plots
#' inp.qc.gt()

inp.qc.gt <- function(inp){

    btween <- ggplot(inp$yg, aes(x=abs(rsnp),y=log(y), group = abs(rsnp), fill = abs(rsnp) )) +
        geom_boxplot(alpha = .2) +
        geom_jitter(width = .05, alpha = .2) +
        guides(fill = "none") +
        theme_bw() +
        labs(
            x = "Genotype rSNP",
            y = "log(Counts)"
        ) +
        ggtitle("Between individual variation")+
        scale_x_continuous(breaks=0:2) +
        scale_y_continuous(limits = c(0, NA))

    y <- data.table(m=rep(inp$gm$m,times=sapply(inp$n,length)),
                        id=rep(1:length(inp$gm$m),times=sapply(inp$n,length)) ,
                        n=unlist(inp$n),
                        p=unlist(inp$p),
                    g=rep(inp$gm$g.ase, times=sapply(inp$n,length))
                    )
    ## add ASE (logit scale ase) n/m or (m-n)/m according if g==-1   
    y[ ,ase:=qlogis(n/m)][g==-1, ase:=-ase]
    dt <- h.inp.qc.gt(y, var="ase", col="no")
    
    if("ai" %in% names(inp)){
        y[,ai:=unlist(inp$ai)][, sdai:=sqrt(unlist(inp$vai))]
        y[,ai.rand:=rnorm(nrow(y),ai,sdai)]
        dt <-rbind(dt, h.inp.qc.gt(y, var="ai", col="yes"))
        sd <- h.inp.qc.gt(y, var="sdai")$V1
        dt[, sd:= c(rep(0,nrow(dt[ref.bias.est=="no",])),sd)]
        ## sample hap estimate from normal dist, when ref.bias.est=="no", sd==0, so no sampling is done
        withCallingHandlers(dt[, rand.est:= rnorm(nrow(dt), V1, sd)],
                            warning=function(w) {
       if (grepl("NAs produced", w$message))
          invokeRestart("muffleWarning")
    } ) 
        
       
        
    }


    
    within <- ggplot(dt, aes(x=as.factor(g), y=rand.est,  color = ref.bias.est)) +
                                        #geom_boxplot(alpha = .2) +
        #geom_point() + 
        geom_jitter(width = .05, shape=1) +
        
        geom_hline(yintercept=0, linetype="dashed") +
        theme_bw() +
        labs(
            x = "Genotype rSNP",
            y = "logit(ASE)"
        ) +
        ggtitle("Within individual variation")

    ## print(length(nb.g))
    
    print(summary(lm(log(y)~abs(rsnp), data=inp$yg))$coefficients)

    print(plot_grid(btween,within))

    return(dt)
}

#' helper to avoid repetition in inp.qc.gt, weigth by prob and merge by individual
#'
#' @param dt
#' @param var
#' @param col creates a new column named "ref.bias.est" with the content in col
#' @keywords helper inp.qc.gt
#' @export
#' @return data table
#' h.inp.qc.gt()

h.inp.qc.gt <- function(dt,var, col=NULL){
    ## weight var by p
    dt2 <- copy(dt)
    dt2[ , eval(var):=get(var)*p]
    
    ## sum var by id
    tmp <- dt2[, sum(get(var)),id]
    tmp <- merge(tmp, unique(dt2[, .(id,g=abs(g))]), by="id")
    ##setnames(tmp, "V1", var)

    if(!is.null(col)){
        tmp[ ,ref.bias.est:=col]
    }
    
    
    return(tmp)
}




#' Extract inputs for a gene-snp association from noGT to make plots for QC purposes
#'
#' @param inp stan input for a gene-snp prepared by stan.trecase.rna.noGT.eff2 (stan.eff.R)
#' @keywords reference panel bias QC
#' @export
#' @return plots
#' inp.qc.nogt()

inp.qc.nogt <- function(inp){


    nb.g <- sapply(inp$NB$p.g,function(j) sum(as.numeric(names(j))*j))
    nb.ngt <- data.table(g=nb.g, counts.log=log(unlist(inp$NB$counts)))

    ## remove 0 counts

    nb.ngt <- nb.ngt[is.finite(rowSums(nb.ngt)),]

    ## print(length(nb.g))
    
    print(summary(lm(g~counts.log, data=nb.ngt))$coefficients)

    betw <- ggplot(nb.ngt, aes(x=g,y=counts.log)) + geom_point(shape=1) +
        geom_smooth(method='lm')+  ggtitle("Between individual variation\n hidden-GT")+ expand_limits(y=0) +theme_bw() +
        labs(
            x = "Genotype rSNP",
            y = "log(Counts)"
        ) 
        

    y <- data.table(m=rep(inp$ase$m,times=sapply(inp$ase$n,length)),
                    id=rep(1:length(inp$ase$m),times=sapply(inp$ase$n,length)) ,
                    n=unlist(inp$ase$n),
                    p=unlist(inp$ase$p),
                    g=unlist(inp$ase$g)
                    )
    ## add ASE (logit scale ase) n/m or (m-n)/m according if g==-1   
    y[ ,ase:=qlogis(n/m)][g==-1, ase:=-ase]
    dt <- h.inp.qc.nogt(y, var="ase", col="no")
    
    if("ai0" %in% names(inp$ase)){
        y[,ai:=unlist(inp$ase$ai0)][, sdai:=sqrt(unlist(inp$ase$vai0))]
        y[,ai.rand:=rnorm(nrow(y),ai,sdai)]
        dt <-rbind(dt, h.inp.qc.nogt(y, var="ai", col="yes"))
        sd <- h.inp.qc.nogt(y, var="sdai")$V1
        dt[, sd:= c(rep(0,nrow(dt[ref.bias.est=="no",])),sd)]
        ## sample hap estimate from normal dist, when ref.bias.est=="no", sd==0, so no sampling is done
        withCallingHandlers(dt[, rand.est:= rnorm(nrow(dt), V1, sd)],
                            warning=function(w) {
       if (grepl("NAs produced", w$message))
          invokeRestart("muffleWarning")
    } )       
        
    }

    
    within <- ggplot(dt, aes(x=as.factor(g.aux), y=rand.est,  color = ref.bias.est)) +
        geom_jitter(width = .05, shape=1) +
        
        geom_hline(yintercept=0, linetype="dashed") +
        theme_bw() +
        labs(
            x = "Genotype rSNP",
            y = "logit(ASE)"
        ) +
        ggtitle("Within individual variation\n hidden-GT")

    print(plot_grid(betw,within))

    return(dt)
}

#' helper to avoid repetition in inp.qc.nogt, weigth by prob and g and merge by individual
#'
#' @param dt, genotypes in column "g"
#' @param var
#' @param col creates a new column named "ref.bias.est" with the content in col
#' @keywords helper inp.qc.nogt
#' @export
#' @return data table
#' h.inp.qc.nogt()

h.inp.qc.nogt <- function(dt,var, col=NULL){
   
    ## weight var by p
    dt2 <- copy(dt)
    dt2[ , eval(var):=get(var)*p]

    ## add aux col for het and hom

    dt2[g==0 | g==2, g.aux:="hom"][abs(g)==1, g.aux:="het"]
    
    ## sum var by  het or hom and id
    tmp <- dt2[, sum(get(var)),.(id, g.aux)]
    ##tmp <- merge(tmp, unique(dt2[, .(id,g=abs(g))]), by="id")
    ##setnames(tmp, "V1", var)

    if(!is.null(col)){
        tmp[ ,ref.bias.est:=col]
    }
    
    
    return(tmp)
}



#' produces table to be added to btrecase.plot
#'
#' @param dt data table with Signif column as it is used in btrecase.plot
#' @param var variables to stratify by, when NULL only Signif
#' @param colors named character vector: values colors and names Signif levels 
#' @keywords table for plot
#' @export
#' @return data table to add to plot, cols are Signif and SNP
#' tab2bplot()

tab2bplot <- function(dt, colors=NULL, var=NULL){
    if(is.null(var)){
        tmp <- dt[, .(SNPs=.N) ,Signif]
    } else {
        tmp <- dt[, .(SNPs=.N) , var]
    }
    
    ## order Signif as is names(colors)
    if(!is.null(colors)){
        tmp <- tmp[order(match(Signif, names(colors)))]
        for(c in names(colors)){
            tmp[Signif==c, color:=colors[c]]
        }
    }
    tmp[,Signif:=factor(Signif, levels(dt$Signif))]
    
    setkey(tmp,Signif)

    return(tmp)
}

    
#' produces table to be added to facet plot for eQTL effect by sign
#'
#' @param dt data table for plotting
#' @param vars character vector with the name of variables to make table
#' @keywords table for plot
#' @export
#' @return data table to add to plot
#' tab2facet()

tab2facet <- function(dt, vars){
    
    g <-  aux.tab2facet(dt, vars)

    tabs <- lapply(1:nrow(g), function(i) {
        
        dt <- dt[get(vars[1])==g[i,Var1] & get(vars[2]) == g[i,Var2] ,.N,.(sign(log2_aFC_mean))]
        dt[, sign:=as.character(sign)][sign=="-1", sign:="<0"][sign=="1", sign:=">0"]
        setkey(dt, sign)
        setnames(dt, c("sign", "N"), c("Effect", "SNPs"))
        return(dt)
    })

    return(tabs)
}


#' Wrapper function to plot Bayesian region of rejections comparing bayesian models
#'
#' @param dt data table for plotting with rej.level column
#' @param colors named character vector: values colors and names Signif levels 
#' @keywords table for plot
#' @export
#' @return data table 
#' aux.tab2facet()

aux.tab2facet <- function(dt, vars){
    tmp <- lapply(vars, function(i) levels(as.factor(dt[[i]])))
    
    g <-  data.table(expand.grid(tmp))

    return(g)
}


#' Format tables for calculating True Discovery Rate
#'
#' @param i table with significant associations (tab2bplot output)
#' @param j name of col with "gold standard"
#' @keywords true discovery rate table
#' @export
#' @return data table with comparison by rejection threshold
#' format.tab

format.tab <- function(i,j) {
    dt <- data.table(Gold_standard = j)
    model <- as.character(i$Signif[which(! i$Signif %in% c("None", "Both", j))])
    if(length(model)) {
        if(length(i[Signif=="Both",SNPs])){
            dt[, eval(model) := i[Signif==model,SNPs] + i[Signif=="Both",SNPs]]
            dt[, Common := i[Signif=="Both",SNPs] ]
            dt[,Total := sum(i$SNPs)]
            dt[, `TDR_%` := round(Common*100/get(model), 2) ]
            dt[,`Power_%` :=round(100*( i[Signif=="Both",SNPs])/(i[Signif==j,SNPs] + i[Signif=="Both",SNPs]), 2)]
            

            if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] + i[Signif=="Both",SNPs]]
            } else {
                dt[, eval(j) :=  i[Signif=="Both",SNPs]]
            }
                      
        
        } else {
            dt[, eval(model) := i[Signif==model,SNPs] ]
            dt[, Common := 0 ]
            dt[,Total := sum(i$SNPs)]
            dt[, `TDR_%` := round(Common*100/get(model), 2) ]
            dt[,`Power_%` := 0]
            
            if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] ]
            } else {
                dt[, eval(j) :=  0]
            }
            
        }
    } else{
        model <- as.character(levels(i$Signif)[which(! levels(i$Signif) %in% c("None", "Both", j))])
        
        if(length(i[Signif=="Both",SNPs])){            
             dt[, eval(model) :=  i[Signif=="Both",SNPs]]
             dt[, Common := i[Signif=="Both",SNPs] ]
             dt[,Total := sum(i$SNPs)]
             dt[, `TDR_%` := round(Common*100/get(model), 2) ]
             dt[,`Power_%` := round(i[Signif=="Both",SNPs]*100/(i[Signif==j,SNPs] + i[Signif=="Both",SNPs]), 2)]
              if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] + i[Signif=="Both",SNPs]]
            } else {
                dt[, eval(j) :=  i[Signif=="Both",SNPs]]
            }
             
        } else {
             dt[, eval(model) :=  0]
             dt[, Common := 0 ]
             dt[,Total := sum(i$SNPs)]
             dt[, `TDR_%` := NA ]
             dt[,`Power_%` := NA]

             if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] ]
            } else {
                dt[, eval(j) :=  0]
            }
        }
           
    }
    
    setcolorder(dt,  c( "Gold_standard", model,  j, "Common", "Total", "TDR_%")) 
    x <- names(i)[which(! names(i) %in% c("Signif", "SNPs", "color"))]
    if(length(x)) {
       
        dt[ , (x):= lapply(x, function(k) unique(i[[k]]))]
       
    }

    dt[, Gold_standard:=NULL]
    return(dt)
}

#' Format tables for calculating sensitiviy and specificity
#'
#' @param i table with significant associations (tab2bplot output)
#' @param j name of col with "gold standard"
#' @keywords true discovery rate table
#' @export
#' @return data table with comparison by rejection threshold
#' format.tab2

format.tab2 <- function(i,j) {
    dt <- data.table(Gold_standard = j)
    model <- as.character(i$Signif[which(! i$Signif %in% c("None", "Both", j))])
    if(length(model)) {
        if(length(i[Signif=="Both",SNPs])){
            dt[, eval(model) := i[Signif==model,SNPs] + i[Signif=="Both",SNPs]]
            dt[, Common := i[Signif=="Both",SNPs] ]
            dt[,Total := sum(i$SNPs)]
            
            dt[,`Power_%` :=round(100*( i[Signif=="Both",SNPs])/(i[Signif==j,SNPs] + i[Signif=="Both",SNPs]), 2)]
            

            if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] + i[Signif=="Both",SNPs]]
            } else {
                dt[, eval(j) :=  i[Signif=="Both",SNPs]]
            }
                      
            
        } else {
            dt[, eval(model) := i[Signif==model,SNPs] ]
            dt[, Common := 0 ]
            dt[,Total := sum(i$SNPs)]
            #dt[, `TDR_%` := round(Common*100/get(model), 2) ]
            dt[,`Power_%` := 0]
            
            if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] ]
            } else {
                dt[, eval(j) :=  0]
            }
            
        }
        
        dt[, `FPR_%` := round(100*( i[Signif==j,SNPs])/(i[Signif==j,SNPs] + i[Signif=="None",SNPs]), 2)]
    } else {
        model <- as.character(levels(i$Signif)[which(! levels(i$Signif) %in% c("None", "Both", j))])
        
        if(length(i[Signif=="Both",SNPs])){            
             dt[, eval(model) :=  i[Signif=="Both",SNPs]]
             dt[, Common := i[Signif=="Both",SNPs] ]
             dt[,Total := sum(i$SNPs)]
             ##dt[, `TDR_%` := round(Common*100/get(model), 2) ]
             dt[,`Power_%` := round(i[Signif=="Both",SNPs]*100/(i[Signif==j,SNPs] + i[Signif=="Both",SNPs]), 2)]
              if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] + i[Signif=="Both",SNPs]]
            } else {
                dt[, eval(j) :=  i[Signif=="Both",SNPs]]
            }
             
        } else {
             dt[, eval(model) :=  0]
             dt[, Common := 0 ]
             dt[,Total := sum(i$SNPs)]
             ##dt[, `TDR_%` := NA ]
             dt[,`Power_%` := NA]

             if(length(i[Signif==j,SNPs])){
                dt[, eval(j) := i[Signif==j,SNPs] ]
            } else {
                dt[, eval(j) :=  0]
            }
        }
        dt[`FPR_%` :=0]
    }
    
    setcolorder(dt,  c( "Gold_standard", model,  j, "Common", "Total", "FPR_%")) 
    x <- names(i)[which(! names(i) %in% c("Signif", "SNPs", "color"))]
    if(length(x)) {
       
        dt[ , (x):= lapply(x, function(k) unique(i[[k]]))]
       
    }

    dt[, Gold_standard:=NULL]
    return(dt)
}


#' Tables of True Disc Rate by Significance and other variables
#'
#' @param dt data table
#' @param var column names with variables to stratify by, including Signif, when NULL only Signif
#' @param cols columns names to stratify TDR table by
#' @param gold name of col with "gold standard"
#' @keywords TRD table
#' @export
#' @return data table  summary
#' TDR.var()

TDR.var <- function(dt, var=NULL, cols, gold){
    tmp <- tab2bplot(dt=dt, var=var)
    ## make character
    tmp[, (cols):=lapply(cols, function(i) as.character(get(i)))]
    comb <-unique(tmp[,(cols),with=F])
    comb <- as.data.table(lapply(comb, function(i) {
        if(is.character(i)){
            i <- paste0("'", i, "'")
        } else {
            i <- i
        }
        
        }))

    tabs <- lapply(1:nrow(comb), function(i){
        txt <- paste0(cols, "==", comb[i,], collapse =" & ")
        tmp[ eval(parse(text=txt)),]
    })

    tabs2 <- rbindlist(lapply(tabs, format.tab, gold))

    return(tabs2)
}

    

#' Recode by categories
#'
#' @param dt data table
#' @param col name of column to discretise
#' @param fun function to apply to col, defaults to NULL
#' @param newcol name of the new column with categories
#' @param cats vector with boundries to categorise "(a,b]" notation
#' @param labs labels for categories
#' @keywords categorise
#' @export
#' @return data table new categorical variable. If value in col out of cats returns NA
#' rec.cat()

rec.cat <- function(dt,col, fun=NULL, newcol,cats,labs) {
    
        for(i in 2:length(cats)){
            if(is.null(fun)){
                dt[get(col) > cats[i-1] & get(col) <= cats[i] , eval(newcol):=labs[i-1]]
            } else {
                txt <- paste0(fun, "(", col, ")")
                dt[eval(parse(text=txt)) > cats[i-1] & eval(parse(text=txt))<= cats[i] , eval(newcol):=labs[i-1]]
            }
        }

    return(dt)
}


#' rasqual vs bayesian plot
#'
#' @param a data table with rasqual output
#' @param b data table with bayesian output by rejection level
#' @param suf character with suffixes to identify data sets ("rasqual", "trec")
#' @param null character vector with col name for null in a and null in b
#' @param rej numeric with rejection levels to consider, if NULL all
#' @keywords compare plot rasqual bayesian
#' @export
#' @return data table with tables and prints btrecase.plot by rasqual fdr for selected rej levels
#'
#' rasq.bay.p()

rasq.bay.p <- function(a,b,suf,null,rej=NULL){
    
    if(is.null(rej)) rej <- unique(b$rej.level)
    dt <- merge(b[rej.level %in% rej,] ,a, by.x=c("Gene_id","tag"), by.y=c("gene_id", "rs_id"), suffixes=paste0(".", suf))
    dt.l <- reshape(dt[,.(Gene_id, tag, rej.level, null.rej, log2_aFC_mean, log2_aFC, null.fdr0.1, null.fdr1, null.fdr5, null.fdr10)],
                    direction="long",
                    varying=list( c("null.fdr0.1", "null.fdr1", "null.fdr5", "null.fdr10")),
                    v.names=null[1],
                    times=as.character(c(0.1, 1, 5, 10)),
                    timevar="Fdr_per")
    names(dt.l)
    if(any(duplicated(null))) null <- paste(null, suf, sep=".")
    dt.l <- add.signif(dt.l, null[1], null[2], suf)
    cols <- c("#999999", "yellow3", "#D55E00", "#0072B2")
    names(cols) <- c("None",suf[1], "Both", suf[2] )
    tab <- tab2bplot(dt.l, colors=cols,
              var=c("Signif", "rej.level","Fdr_per"))
    names(dt.l)
    for(k in rej){
        tmp <- dt.l[rej.level == k,]
        
        tab.rej <- tab[rej.level == k,]
        
        lev <- unique(tab.rej$Fdr_per)
        tmp[ , Fdr_per:=factor(Fdr_per, levels=lev)]
        
        tab.fdr <- lapply(lev, function(i)  tab.rej[Fdr_per == i,])
        
        gl <- lapply(tab.fdr, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = 14,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))), xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))

        dt.tab <- data.table(Fdr_per=factor(lev), grob=gl)

        lab <- paste0('FDR = ', lev, '%')
        names(lab) <- lev

        

##' # Rasqual with different FDR relative to Trec
#+ fig.width= 12, fig.height=21
        p <- btrecase.plot(dt=tmp[Signif != "None" ,] , x2=c(rep("log2_aFC_mean",3),null[2]) ,
                           x1=c(rep("log2_aFC",3), null[1]),
                           yl=paste("eQTL effect ", suf[2]),
                           xl=paste("eQTL effect ", suf[1]),
                           col=suf,
                           title=paste(paste(suf, collapse=" vs "), "\n for rejection zone = ", k),
                           title.size=16,
                           axis.title = 14,
                           axis.text=12,
                           legend.title=14,
                           legend.text=12,
                           legend.symbol=5,
                           point.size=3
                           ) +  facet_grid(Fdr_per~., labeller=labeller(Fdr_per=lab)) +
            theme(strip.background=element_rect(fill="white") ,strip.text.y = element_text(size = 14)) +
            geom_custom(data=dt.tab, aes(grob=grob), x = 0.15, y = 0.75)
        print(p)
    }
    
    return(tab)

 
}


#' qqplot for mixture of normals
#'
#' @param data, vector with observed data
#' @param normals list with 3 elements: mu, mean of each gaussian, sd, sd of each gaussian, mix, mixing proportion
#' @param p vector with probabilities to get quantiles, defaul to 10^4 values
#' @param plot, whether to return plot or data table to plot, defaults to plot
#' @keywords mixture qqplot gaussians
#' @export
#'
#' qq.mix()
#' 

qq.mix <- function(data,normals,p=seq(0,1,0.0001), plot=NULL){
    x <- rmvnormmix(n=length(data), lambda=normals$mix, mu=normals$mu, sigma=normals$sd)

    dt.q <- data.table(dat = quantile(data, probs=p), x= quantile(x, probs=p ))
    if(!is.null(plot)) return(dt.q)
    
    p <- ggplot(dt.q, aes(x,dat)) +
        geom_point() +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="red") + ggtitle(paste("Mixture of", unique(sapply(normals,length)),  "Gaussians")) +
        ylab("Observed quantiles") +
        xlab("Mix-normal quantiles")

    return(p)
}


#' qqplot for comparing priors
#'
#' @param p1 prior 1: list with 3 elements: mu, mean of each gaussian, sd, sd of each gaussian, mix, mixing proportion
#' @param p2 list with 3 elements: mu, mean of each gaussian, sd, sd of each gaussian, mix, mixing proportion
#' @param n number of times to sample from p1 or p2, defaulst to 1000
#' @param p vector with probabilities to get quantiles, defaul to 10^4 values
#' @keywords mixture qqplot gaussians prior
#' @export
#'
#' qq.p()

qq.p <- function(p1,p2,p=seq(0,1,0.0001), n=1000){


    com.p <- lapply(list(p1,p2), function(i) {
        quantile(as.numeric(rmvnormmix(n=n, lambda=i$mix, mu=i$mu, sigma=i$sd)),
                 probs=p)
    })
        
    comp.p <- Reduce(function(a,b) data.table(px=a, py=b), com.p)

    
    p <- ggplot(comp.p, aes(px,py)) +
        geom_point() +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="red") + ggtitle(paste("Mixture of", unique(sapply(p1,length)),  "Gaussians\n spike with sd = 0 or sd =", p2$sd[1])) +
        ylab(paste0("Prior with sd=",p2$sd[1])) +
        xlab("Prior with sd=0")

    return(p)
}


#' plot and tables comparing freq and bayesian methods, wrapper for tab2bplot and btrecase.plot
#'
#' @param a list with each element a merge of freq and bayesian
#' @param fac named list with each element the values for each factor to split summary, and name the column name for the factor in a. Defaults to Signif, input for tab2bplot. All elements of fac should have the same length
#' @param colors named character vector: values colors and names Signif levels
#' @param var named characeter vector:name is the col name to subset each element of a, value is the value of var to subset by,defaults to NULL
#' @param title character, names of a will be added at the end for reference
#' @param xpos x position to add table into graph
#' @param ypos y position to add table to plot
#' @param facet.fac one or 2 factors to facet ggplot, defaults to fac names
#' @param x1 names of cols with mean, CI min and max and null for condition1
#' @param x2 names of cols with mean, CI min and max and null for condition2
#' @param s number of row to sample when data table is too long, defaults to 5000
#' @param rx vector with range for x axis
#' @param ry vector with range for y axis
#' @param xl character with xlab
#' @param yl character with ylab
#' @param col when usimg color for symbol based on sig or not sig, provide identifier for first and second column being significant
#' @param title.size size for title
#' @param axis.title size for axis.title 
#' @param axis.text size for axis.text.x and axis.text.y 
#' @param legend.title size for legend title
#' @param legend.text size for legend.text
#' @param legend.symbol size for legend symbols
#' @param tab.size size for text in table
#' @param point size for geom_points
#' @param plot prints plot by default unless not null
#' @export
#' @return prints plot and return tables for each element of "a" by "fac"
#' plot_tab()

plot_tab <- function(a,fac=NULL, colors, var=NULL, title=NULL, xpos, ypos, facet.fac=NULL,x1, x2, s=5000, rx=NULL,ry=NULL,xl, yl, col=NULL, title.size=12,axis.title=12, axis.text=10, legend.title=12, legend.text=10, legend.symbol=4, tab.size=16, point.size=3, plot= NULL){
    
    l <-  mapply(function(a, b) {
        table <- tab2bplot(dt=a, var= c("Signif", names(fac)), colors=colors)
       
        tables <- lapply(1:length(fac[[1]]), function(i) {
            expl <- lapply(names(fac), function(j) paste(j, "==", fac[[j]][i]))
            exp <- Reduce(function(p,q) paste(p,q, sep=" & "), expl)
            dt <- table[eval(parse(text=exp)),]
            return(dt)
        })

        
        
        gl <- lapply(tables, function(i) tableGrob(i[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = tab.size,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(i$color, rep("black", 4)))))))#, xmin=-2.5, xmax=-1.5, ymin=1.5, ymax=2.5)))


        dt.tab <- data.table(grob=gl)

        if(is.null(facet.fac)) {
            dt.tab[,names(fac):=fac ]
        } else {
            dt.tab[, eval(facet.fac):=fac[facet.fac]]
        }
        
        if(is.null(plot)) {

            tit <- if(!is.null(title)) paste(title, b) else NULL

            if(!is.null(var)) a <- a[get(names(var)) == var,]

            
            p <- btrecase.plot(dt=a[Signif != "None" ,] , x1=x1 ,
                               x2=x2,
                               s=s,
                               xl=xl,
                               yl=yl,
                               col=col,
                               title=tit,
                               title.size=title.size,
                               axis.title = axis.title,
                               axis.text=axis.text,
                               legend.title=legend.title,
                               legend.text=legend.text,
                               legend.symbol=legend.symbol,
                               point.size=point.size
                               )

            if(is.null(facet.fac)) facet.fac=names(fac)

            if(length(facet.fac) == 1 ){
                exp <- paste(facet.fac, "~.")

            } else {
                exp <- paste(facet.fac, collapse=" ~ ")
            }
            

            p <- p + facet_grid(as.formula(exp)) + #, labeller=labeller(rej.level=lab))+
                theme(strip.background=element_rect(fill="white") ,strip.text.y = element_text(size = 14), strip.text.x = element_text(size = 14)) +
                geom_custom(data=dt.tab, aes(grob=grob), x = xpos, y = ypos)
            

            
            print(p)
        }
        
        return(tables)
    } ,
    a=a,
    b=names(a),
    SIMPLIFY=F)

    return(l)
    
}
 
#' tables comparing freq and bayesian methods, wrapper for tab2bplot
#'
#' @param a list with each element a merge of freq and bayesian
#' @param fac named list with each element the values for each factor to split summary, and name the column name for the factor in a. Defaults to Signif, input for tab2bplot. All elements of fac should have the same length
#' @export
#' @return tables for each element of "a" by "fac"
#' only_tab()

only_tab <- function(a,fac=NULL){
    
    l <-  mapply(function(a, b) {
        table <- tab2bplot(dt=a, var= c("Signif", names(fac)), colors=NULL)

        tables <- lapply(1:length(fac[[1]]), function(i) {
            expl <- lapply(names(fac), function(j) paste(j, "==", fac[[j]][i]))
            exp <- Reduce(function(p,q) paste(p,q, sep=" & "), expl)
            dt <- table[eval(parse(text=exp)),]
            return(dt)
        })
        
       
        return(tables)
    } ,
    a=a,
    b=names(a),
    SIMPLIFY=F)

    return(l)
    
 }

#' aux for plotting bayesian models by null
#'
#' @param a stan sum for model1
#' @param b stan sum model2
#' @param null name of the col with null in both datasets
#' @param byx cols in a to merge a and b
#' @param byy cols in b to merge a and b
#' @param suffixes for merging
#' @param siglevel significant level, defaults to 95
#' @param colsig character with the names to call significant hits in each dataset
#' @param xl title for x axes
#' @param yl title for y axes
#' @param title graph title
#' @param size.tab text size for table in graph
#' @param xmin min value for location of table in graph
#' @param xmax max value for location of table in graph
#' @param ymin min value for location of table in graph
#' @param ymax max value for location of table in graph
#' @param prior name of prior used in the model to include in title, deafaults to NULL
#' @export
#' @keywords plot comparing model
#' @return plot made by btrecase.plot
#' merge.plot()
 
merge.plot <- function(a,b, null="null.95", byx, byy, suffixes,siglevel=95, colsig, xl, yl, title, size.tab, xmin, xmax, ymin, ymax, prior=NULL){
    dt <- merge.aux(a,b,null, byx, byy, suffixes, colsig)
    
    cols <- setNames( c("#999999", "yellow3","#0072B2",  "#D55E00") , c("None", colsig, "Both"))
    tab <- tab2bplot(dt, colors=cols)
    ## get CI lims
    cilims <- 50+as.numeric(siglevel)/2*c(-1,1)
    if(!is.null(prior)) title <- paste0(title, ", prior = ", prior)
    x=c("log2_aFC_mean", paste0("log2_aFC_", cilims, "%"), null)
    p <- btrecase.plot(dt=dt[Signif !="None",],
                       x1=paste0(x, suffixes[1]),
                       x2=paste0(x, suffixes[2]),
                       xl=xl, yl=yl, col=colsig, title=title)
    p <- p +  annotation_custom(tableGrob(tab[,.(Signif,SNPs)], rows=NULL, theme=ttheme_minimal(base_size = size.tab,padding = unit(c(2, 1.5), "mm"), core=list(fg_params=list(col=c(tab$color, rep("black", 4)))))), xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    return(p)
}

#' aux for merge.plot
#' @param a stan sum for model1
#' @param b stan sum model2
#' @param null name of the col with null in both datasets
#' @param byx cols in a to merge a and b
#' @param byy cols in b to merge a and b
#' @param suffixes for merging
#' @param colsig character with the names to call significant hits in each dataset
#' @export
#' @return data table 
#' merge.aux()

merge.aux <- function(a,b,null, byx, byy, suffixes, colsig){
    dt <- merge(a,b,by.x=byx, by.y=byy, suffixes=suffixes)
    dt <- add.signif(dt, paste0(null, suffixes[1]), paste0(null, suffixes[2]), col=colsig)
    return(dt)
}



#' helper function to calculate true positives and power by FDR for Rasqual, DEseq and Btrecase (GT) for comparing to whole GEUVADIS
#' @param total data table with total sig calls (prepared in geu_Figs.R)
#' @param common data table with common sig calls (prepared in geu_Figs.R)
#' @param fdr data table with btrecase fdr, columns post.level and FDR, other cols allowed
#' @param N number of positive associations in GEU (for the associations tested in the 3 methods)
#' @export
#' @return data table with relevant cols
#' help.geu()

help.geu  <- function(total, common, fdr, N){
    comp.geu <- merge(common, total, by=c("Method", "Fdr", "PIP"), suffixes=c(".common", ".total"))
    comp.geu <- merge(comp.geu, fdr[, .(post.level,FDR)], by.x="PIP", by.y="post.level", all.x=T)
    comp.geu[is.na(Fdr), Fdr:=FDR][, FDR:=NULL][ , TP:= round(N.common/N.total, 2)]

    ## add total number of postive associations  in GEU (same for all methods and FDRs)

    comp.geu[, sig.geu:=N ]
    comp.geu[ , `Power_%`:=round(100*N.common/sig.geu,2)]
    comp.geu[, `TDR_%`:= TP*100]
    comp.geu[ ,Discoveries:=N.total]
    return(comp.geu)
}

#' add r2 to stan summary to each rSNP relative to top rSNP
#' 
#' @param sum summary with stan output
#' @param gene gene_id to add r2
#' @param top name of column with info for top rSNP,
#' @param notop character with value in top column for SNPs that are not the top SNP, defaults to 0
#' @param r2 matrix with colnames rSNPs and rownames rSNPs, previoulsy computed matrix with r2 for all combination of rSNPs. rSNPs run for gene.
#' @param col name for column with r2 value, defaults to r2
#' @export
#' @return data table for gene only (subset of stan summary for gene=gene) with r2 col
#' r2.stan()

r2.stan  <- function(sum, top="Top", notop="0", gene , r2, col="r2"){
    dt <- sum[Gene_id ==gene,]
    top.r <- dt[get(top) != notop, tag ]
    rest.r <- dt[get(top) == notop, tag ]
    r2.sub <- as.data.table(r2[c(top.r, rest.r),top.r, drop=F],keep.rownames=T)
    setnames(r2.sub, names(r2.sub), c("tag", col))
    dt <- merge(dt, r2.sub, by="tag", sort=F)
    return(dt)
}

    
#' get rsid for a SNP
#'
#' @param rsid file with variant information formatted as ensembl ftp.ensembl.org/pub/grch37/current/variation/gvf/homo_sapiens/homo_sapiens-chr22.gvf.gz for appropiate built and chromosome, defaults to none
#' @param snps vector with POS:REF:ALT to look snp for. POS must be in same built as rsid. Accepts many snps but for same chrom (rsid file)
#' @export
#' @keywords rsid SNP
#' @return data table with POS:REF:ALT and rsid
#' get.rsid

get.rsid <- function(rsid, snps){
    pos <- as.numeric(sub(":[A-Z].*", "\\1", snps))
    ref <- sub(".*:(.+):.*", "\\1", snps)
    alt <- sub(".*:(.+)", "\\1", snps)
    ## make grep -e syntax allowing multiple pos
    ex <- paste(paste("-e", pos), collapse=" ")
    dt <- fread(cmd=paste("zgrep", ex, rsid))

    ## select data from dt

    m <- mapply(function(a,b,c) {
        l <- list(dt[a == dt$V4, which=T],
                  grep(paste("Reference_seq", b, sep="="), dt$V9),
                  c(grep(paste0("Variant_seq=",c , "[;,]"), dt$V9), grep(paste0("Variant_seq=.*,",c, ";"), dt$V9))
                  )
        r <- Reduce(intersect, l)
        if(length(r) ==1){
            rs <- dt[r, sub(".*Dbxref=dbSNP_151:(.+);evidence.*", "\\1", V9)]
            return(data.table(SNP=paste(a,b,c, sep=":"), rsid=rs))
        } else {
            return(data.table(SNP=paste(a,b,c, sep=":"), rsid="-"))
        }
    },  
    a=pos,
    b=ref,
    c=alt,
    SIMPLIFY=F)
    return(rbindlist(m))
}

#' convert to similar fdr for plotting
#'
#' @param dt data table with Method and Fdr col
#' @param fdr vector with new values for fdr to keep consistency across groups
#' @param method name of method to change fdr values, defaults to Btrecase
#' @export
#' @return data table with new col Fdr.approx
#' approx.fdr()
#'
approx.fdr <- function(dt, fdr=c(0.001, 0.01, 0.05, 0.1), method="Btrecase"){
    dt <- dt[ Method==method | (Method!=method & Fdr %in% fdr) ,]
    setkey(dt, Method, Fdr)
    return(dt[,Fdr.approx:= rep(fdr, length(unique(dt$Method)))])
}


#' fdr plot for exernal validity, using some libraries/commands from Chris
#' 
#' @param dt data table with summaries and relevant cols
#' @param labels for legend, defaults to NULL
#' @param x name for x col, defaults to TDR_%
#' @param y name for y col, defaults to Power_%
#' @param lab.col, name of col to add labels, defaults to "Discoveries"
#' @param type character specifying if plotting by associations or genes
#' @param col name of col to get gold standard associations or egenes
#' @param path whether to link observations by method, defaults to yes
#' @param label whether to add labels, defaults to yes
#' @export
#' @return ggplot object
#' fdr.plot()
#'
fdr.plot <- function(dt, x="TDR_%", y="Power_%", lab.col="Discoveries",  labels=NULL, type=c("associations", "eGenes"), col, path="yes", label="yes"){
    ##aux function for consistent scales
    scaleFUN <- function(x) sprintf("%.2f", x)  ## for making consistent scales across plots
    if(x== "TDR_%") { ## make compatible with cols already in right scale
        dt[, tdr:=`TDR_%`/100][, power:=`Power_%`/100]
        x="tdr"
        y="power"
    }
    
    p <- ggplot(data=dt, aes(get(x),  get(y))) +
        lemon::geom_pointpath(aes(get(x),  get(y),col=Method,group=Method,pch=Method),lwd=2)+
        theme_cowplot(font_size=10)
    if(path=="yes"){
        p <- p +
        ##geom_path(mapping= aes(`TDR_%`/100,  `Power_%`/100, color=Method), linetype="dotted") +
            geom_path( mapping=aes(get(x),  get(y), group=as.factor(Fdr.approx)),color="gray", linetype="dashed",lwd=0.5)
    }
    if(label=="yes"){
        p <- p +
        ##geom_point(mapping= aes(`TDR_%`/100,  `Power_%`/100, color=Method, shape=Method)) +
            geom_label_repel(mapping= aes(get(x),  get(y), label=get(lab.col), color=Method ),show.legend=FALSE,size=2 )
    }
    p <- p +
        xlab("Positive predicted value") +
        ylab("Sensitivity") +
        ggtitle(paste0('"True" ',type, ': ', unique(dt[[col]]))) +      
        scale_y_continuous(labels=scaleFUN)
    if(is.null(labels)){
        p <- p + scale_colour_ipsum()
    } else {
        
        p <- p +  scale_colour_ipsum(labels=lab) + scale_shape_discrete(labels=lab)
    }
    return(p)

}


## code for plot with 2 legends at the bottom
## fdr.pa <- lapply(fdr.tabs, function(i) {
##     ggplot(data=i) + #, aes(`TDR_%`/100,  `Power_%`/100, color=Method, shape=Method, label=Discoveries )) +
##         geom_path(mapping= aes(`TDR_%`/100,  `Power_%`/100, color=Method), linetype="dotted") +
##         geom_path( mapping=aes(`TDR_%`/100,  `Power_%`/100, linetype=as.factor(Fdr.approx)) ,color="gray") +
##         geom_point(mapping= aes(`TDR_%`/100,  `Power_%`/100, color=Method, shape=Method)) +
##         geom_text_repel(mapping= aes(`TDR_%`/100,  `Power_%`/100, label=Discoveries, color=Method ) ) +
##         xlab("Positive predicted value") +
##         ylab("Sensitivity") +
##         ggtitle(paste0('"True" associations: ', unique(i[["Gtex-ebv"]]))) +
##         ## change label in legend
##         scale_colour_manual(values=c("#009E73", "#CC79A7", "blue"), labels=lab) +
##         scale_shape_discrete(labels=lab) +
##         scale_y_continuous(labels=scaleFUN)  +
##         scale_linetype_manual(name="FDR",
##                               values=c("solid","twodash", "dotted",  "dashed"))  +
##         theme(legend.position="bottom", legend.box="vertical", legend.spacing.y=unit(0, 'cm')) +
##         guides(color = guide_legend(order=1),
##                linetype = guide_legend(order=2),
##                shape=FALSE)          #xlim(0.15,1)
## })

