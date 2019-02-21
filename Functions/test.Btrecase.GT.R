library(testthat)
source("/home/ev250/Bayesian_inf/trecase/Functions/Btrecase.GT.R")



test_that("input files except covariates", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=5*10^5,
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     model="both")
  
  expect_that(inp, throws_error())
  
})

test_that("numeric inputs", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=5*10^5,
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets="5",min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/")
  
  expect_that(inp, throws_error("invalid arguments:"))
  
})

test_that("numeric snps and model='both'", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=0,
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                      
                     prefix="num.snps.both.test")
  
  expect_that(inp, is_a("data.table"))
  
})

test_that("numeric snps and model='trec'", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=0,
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     prefix="num.snps.trec.test",
                     
                     model="trec")
  
  expect_that(inp, is_a("data.table"))
  
})

test_that("numeric snps and model='trecase'", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=0,
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     prefix="num.snps.trecase.test",
                     
                     model="trecase")
  
  expect_that(inp, is_a("data.table"))
  
})

test_that("selected snps and model='both'", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=c("50315933:G:C", "50316471:A:G","50312651:C:T", "50312672:A:G"),
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold=.9,q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     prefix="char.snps.both.test")
  
  expect_that(inp, is_a("data.table"))
  
})


test_that("selected snps, no tagging and model='both'", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=c("50315933:G:C", "50316471:A:G","50312651:C:T", "50312672:A:G"),
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold="no",q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     prefix="char.snps.notag.both.test")
  
  expect_that(inp, is_a("data.table"))
  
})

test_that("numeric snps, tagging, quick test and model='both'", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps=0,
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold=0.9,q.test="yes",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     prefix="num.snps.qtest.both.test")
  
  expect_that(inp, is_a("data.table"))
  
})

test_that("selected snps,no tagging, snp not in reference panel, and model='both'", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps='50315933:G:C',
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold="no",q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     prefix="char.snps.not.in ref.pan.both.test")
  
  expect_that(inp, is_a("data.table"))
  
})

test_that("selected snps,no tagging, hap not in reference panel, and model='trecase' ", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps='50318735:C:T',
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold="no",q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     model='trecase',
                     
                     prefix="char.snps.hap.not.in ref.pan.trecase.test")
  
  expect_that(inp, is_a("character"))
  
})

test_that("selected snps,no tagging, hap not in reference panel, and model='trecase' ", {
  
  inp <- btrecase.gt(gene="ENSG00000184164",
                     chr=22,
                     snps='50312404:C:T',
                     counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.GTqc.txt',
                     
                     covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.GT.rds',
                     
                     e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt',
                     
                     gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt',
                     
                     vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.GTqc.vcf.gz' ,
                     
                     le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz' ,
                     
                     h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz',
                     
                     nhets=5,min.ase=5,min.ase.het=5,tag.threshold="no",q.test="no",
                     
                     out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/",
                     
                     model='trecase',
                     
                     prefix="char.snps.hom.trecase.test")
  
  expect_that(inp, throws_error("Missing GT or homozygous snps in all samples"))
  
})



ex=fread("/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/num.snps.both.test.trecase.excluded.rsnps.txt")


gene="ENSG00000100243"
chr=22
snps=5*10^5
counts.f='/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/b37_filtered.raw_counts.txt'

covariates= '/mrc-bsu/scratch/ev250/EGEUV1/quant/RNA_counts/library.size.comp.rds'

e.snps= '/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/chr22.fSNPs.txt'

gene.coord='/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input/gene_data_longest.exons.txt'

vcf='/mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/chr22.ASE.allsamples.vcf.gz'

le.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz'

h.file='/home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz'

nhets=5
min.ase=5
min.ase.het=5
tag.threshold=.9
q.test="no"
population="EUR"
prob=NULL
model="trecase"
## alternative

prob=0.9
tag.threshold="no"

## out="/mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/output/chr22/"
# 
# prefix="num.snps.both.test"
# 
# model="both"


################## emedlab #############
gene="ENSG00000000419"
chr=20
snps=5*10^5
counts.f='/mrc-bsu/scratch/ev250/emedlab/RNA/RNA_counts/groups/Genentech_Paired_Baseline_Synovium_RA_Geno.txt'


vcf='/mrc-bsu/scratch/ev250/emedlab/DNA/imputed/PEAC_1000G_Phase3_Oct14_chr20.vcf.gz'


nhets=5

tag.threshold=.9
