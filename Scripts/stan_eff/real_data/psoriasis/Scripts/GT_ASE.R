library(data.table)

#' Function to prepare body for vcf merging GT and ASE fields
#'
#' Format individual samples per chr with GT:AS  as required for vcf
#' @param GT full name to file with GT vcf body 
#' @param ASE full name to file with ASE counts, ouput from rule ASE
#' @param out full name for output file 
#' @keywords vcf GT ASE
#' @export
#' @return saves a tab delimited with vcf body per sample per chr 
#' GTaseVcf()

GTaseVcf <- function(GT, ASE, out){
    gt <- fread(GT)
    ase <- fread(ASE)
    dt <- merge(gt, ase[,c("contig","position","refAllele","altAllele","refCount","altCount","totalCount"),with=F] ,
                by.x=paste0("V",c(1:2,4:5)),
                by.y=c("contig","position","refAllele","altAllele"), all.x=T)
    ## replace NA with 0
    for (i in names(dt)){
        dt[is.na(get(i)), (i):=0]
    }
    ## if V9 is 1|0, swap refCount with altCount ##
    dt[V9=="1|0", refCount:=altCount]
    dt[V9=="1|0", altCount:=totalCount-altCount]

    ## merge refCount with altCount
    dt[,V8:=rep("GT:AS",nrow(dt))]
    dt[,V9:=paste0(V9,":",refCount,",",altCount)]

    ## add info col
    dt[, INFO:="."]

    ## remove columns and reorder
    dt[,c('refCount', 'altCount', 'totalCount'):=NULL]
    setcolorder(dt, c(names(gt)[1:7], "INFO", names(gt)[8:9]))
    setkey(dt, V2)
     
    ## save as tab del
    write.table(dt, out,row.names=F,col.names=F,quote=F,sep="\t")
     
}

GTaseVcf(snakemake@input[['GT']],snakemake@input[['ASE']], snakemake@output[[1]])
