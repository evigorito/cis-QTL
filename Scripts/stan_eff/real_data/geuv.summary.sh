########################################################
## Apply bayesian trecase to GEUV data of EUR ancestry
########################################################

## module load samtools-1.4-gcc-5.4.0-derfxbk

###### functions ######
. /home/ev250/Cincinatti/Functions/various.sh

###################################################
## select samples to download based on EU ancestry

cd /mrc-bsu/scratch/ev250/EGEUV1/sample_info
wget -k  ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/GEUV/E-GEUV-1/E-GEUV-1.sdrf.txt

## From file select  samples with GBR ancestry (EUR) and get ftp address for fasq files

## I had bad experience with bam files in the past, was adviced to start from fasq files.

## get col number for ftp
col=$(awk -v RS='\t' '/Comment\[FASTQ_URI\]/{print NR; exit}' E-GEUV-1.sdrf.txt)

grep GBR E-GEUV-1.sdrf.txt | cut -f 1,$col  > GBR.samples

## some samples are already downloaded in /mrc-bsu/scratch/ev250/EGEUV1/RNA_seq_fastq/
## use samples names of available samples to list samples to download
cat /scratch/ev250/EGEUV1/RNA_seq_fastq/sample_names.txt |cut -d ' ' -f1 > samples.in
grep -v -f samples.in GBR.samples > samples.to.download


## download RNA-seq samples
cd /mrc-bsu/scratch/ev250/EGEUV1/RNA_seq_fastq
cat /mrc-bsu/scratch/ev250/EGEUV1/sample_info/samples.to.download |cut -f2 | parallel --gnu "wget {}"

## download DNA genotype data: start with chr22
cd /mrc-bsu/scratch/ev250/EGEUV1/DNA

wget -r --no-parent -A"GEUVADIS.chr22.*" -nH --cut-dirs=6 http://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/genotypes/


######################################################################
## RNA-seq: align fasq files using STAR in build37

## some samples are already aligned and stored in /mrc-bsu/scratch/ev250/EGEUV1/quant/STAR/built37, I need to align the ones I downloaded (/mrc-bsu/scratch/ev250/EGEUV1/sample_info/samples.to.download)

## I adapt the code I run before for aligment /home/ev250/Genotyping_RNA_seq/Scripts/star_map_sub.sh and star_map.sh and save it in '/home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/star37'


########################################################################
## DNA: extract relevant samples from file 

cd /mrc-bsu/scratch/ev250/EGEUV1/DNA

## issues with DNA chr22 header, as I got before with chr14. Need to add:

##contig=<ID=22,assembly=b37,length=51304566>
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihoods">
##FORMAT=<ID=PP,Number=1,Type=Float,Description="PP">
##FORMAT=<ID=BD,Number=1,Type=Float,Description="BD">

## length chr22
https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
chr22	51304566

## modify header (sed '$i inserts line before the last)

bcftools view -h GEUVADIS.chr22.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz | sed -e '3i\##contig=<ID=22,assembly=b37,length=51304566>' | sed 's/FORMAT=<ID=GL,Number=./FORMAT=<ID=GL,Number=G/' | sed '$i\##FORMAT=<ID=PP,Number=1,Type=Float,Description="PP"> ' | sed '$i\##FORMAT=<ID=BD,Number=1,Type=Float,Description="BD"> '  >  GEUVADIS.chr22.mod.header.vcf

## append vcf body to new header
bcftools view -H GEUVADIS.chr22.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz >> GEUVADIS.chr22.mod.header.vcf

## extract relevant samples
cat /mrc-bsu/scratch/ev250/EGEUV1/sample_info/GBR.samples | cut -f1 | uniq > /mrc-bsu/scratch/ev250/EGEUV1/sample_info/GBR.sample.id

bcftools view -S  /mrc-bsu/scratch/ev250/EGEUV1/sample_info/GBR.sample.id > GEUVADIS.chr22.GBR.vcf

## vcf is already phased, dont need to run shapeit

## compress and index

bgzip -c GEUVADIS.chr22.GBR.vcf > GEUVADIS.chr22.GBR.vcf.gz
tabix -p vcf GEUVADIS.chr22.GBR.vcf.gz


###########################################################################################
## Use RNA bam files and DNA shapeit file for ASE

## Run phaser.array.sh and GBR22.sh, each job one sample, chr22.

########################################################################################
## Prepare inputs for eQTL analysis
####################################

################################################################
## extract GT information from each sample and merge it with ASE information

vcf4AS /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE 22 22

## format in inputs.R: save samples in "...for.AS.tab"
## convert tab files into vcf
tab2vcf /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE for.AS.tab

## merge all samples into 1 vcf per chr

mergevcf /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE 22 22

## rm intermediate files: not done


## QC vcf file, make sure GT is in the right format: inputs.R
## make sure counts and vcf have the same samples:inputs.R

################################### Extract snp coordinates with REF and ALT alleles

vcf4rasqual /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE 22 22 ASE.allsamples.vcf.gz  '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' /mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input snp.coord.txt

################# Prepare counts per gene and format inputs for Btrecase #######

## done in /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/QTL_input/inputs.R

################## Run model for chr22 ##########

## Started in 'top.snp.set.22.R' for gene-snp associations based on Chris eqtl geuvadis data.
## In this script I call others to run the top snp but also a snp in low LD to set rules to avoid running unecessary tests.

## Expanding into chr22.R, in inputs.R I could make a list of genes to test based on chr and level of expression so I can then make an array variable that goes through that list to run eQTL.

#########################################################################################
################## Discrepancy between phasing on ref panel and GEUVADIS DNA-download ###
#########################################################################################

## We were recommended to use 1000GP phase3 so I now select samples from ref panel.

## /home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3 has the following files for chr22:
chr22.bcf.gz  chr22.bcf.gz.csi  chr22.hap.gz  chr22.legend.gz  chr22.samples

## chr22.bcf.gz doesnt have alt allele, open the file in R, add it and re-make bcf, for some reason the convertion from hap/sam/legend files does not add some of the basic fields. 

## Get samples of interest: some arent in vcf, --force-samples exclude them w/o error message

## get header but remove everything except GT
bcftools view  -S /mrc-bsu/scratch/ev250/EGEUV1/sample_info/GBR.sample.id /home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/chr22.bcf.gz --force-samples -h > /mrc-bsu/scratch/ev250/EGEUV1/DNA/chr22.1KGP_P3.vcf

## get body with GT to open in R
bcftools view -S /mrc-bsu/scratch/ev250/EGEUV1/sample_info/GBR.sample.id /home/ev250/rds/rds-cew54-wallace-share/Data/reference/1000GP_Phase3/chr22.bcf.gz --force-samples -Ou | bcftools query -f "%CHROM %POS %ID %REF %ALT %QUAL %FILTER %AC %AN[ %GT]\n" > /mrc-bsu/scratch/ev250/EGEUV1/DNA/chr22.1KGP_P3.txt

## fixed ALT allele and select SNPs only, also phaser requires FILTER column in vcf to be PASS, changes in /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/QTL_input/fix.alt.R

## append vcf body to new header
cat /mrc-bsu/scratch/ev250/EGEUV1/DNA/chr22.1KGP_P3body.txt >> /mrc-bsu/scratch/ev250/EGEUV1/DNA/chr22.1KGP_P3.vcf

## compress and index
bgzip -c /mrc-bsu/scratch/ev250/EGEUV1/DNA/chr22.1KGP_P3.vcf > /mrc-bsu/scratch/ev250/EGEUV1/DNA/chr22.1KGP_P3.vcf.gz
tabix -p vcf /mrc-bsu/scratch/ev250/EGEUV1/DNA/chr22.1KGP_P3.vcf.gz

## Run phaser.array.sh and GBR22.sh, each job one sample, chr22
## re run pipeline
##

## I compared lm. gt and nogt for chr22. We changed the prior for noGT to be bj~N(0,0.2) and re-run. All seems ok, analysis in comp.out.chr22.GT.noGT.lm.R

## Next steps:

##1) run noGT with genotypes called from RNA

##2) improve input preparation, need to use bedtools to count only once reads overlapping more than one snp.


############################################################################################
## 1) ##### Genotyping from RNA #########
############################################################################################

## start with chr22: tried all chrs but calling variants is too slow. Better by chr.
## run in call_var_rna.chr22.sh: /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.vcf, calls variants, changed samples names and added genomic annotations.

## extract FORMAT/INFO fields for chr 22 to look at variant calling in R

bcftools query -f  '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%ANN[\t%GT\t%DP\t%AD\t%RPB\t%MQ\t%VDB\t%MQB\t%BQB]\n' /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.vcf > /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.txt

bcftools query -f  '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%ANN[\t%GT\t%DP\t%AD\t%RPB\t%MQ\t%VDB\t%MQB\t%BQB]\n' /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.vcf -H | head -1 > /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.header.txt

## compare RNA genotyping with DNA genotyping in /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/call_var/RNA.genotyping.QC.R:

# 1) Select variants called in DNA and use DP per sample >=10

# 2) Prepare vcf for shapeit.

## header with no INFO and FORMAT GT only,  exclude samples not in DNA (RNA.genotyping.QC.R)

## needed to do this order of commands: -s didnt work in annotate and bcftools view adds AC and AN in info and then filled not sure how when added vcf body. To avoid it use -I

bcftools view -s ^HG00104,HG00124,HG00134,HG00135,HG00152,HG00156,HG00247,HG00249 /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20.vcf -h -I | bcftools annotate  -x INFO,^FORMAT/GT > /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20_GTonly.vcf

## body saved in /RNA.genotyping.QC.R

## append vcf body to new header
cat /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20_GTonly.txt >>  /mrc-bsu/scratch/ev250/EGEUV1/quant/call_vars/b37/HG96_2215_chr22_q20_GTonly.vcf


## shapeit: run in rna_chr22.sh

## need to exclude missing values as I get error (more than 10% of variants missing in any sample). Try phasing each sample after removing missing values. I tried:

## bcftools view -s HG00096 -e 'GT="./."' HG96_2215_chr22_q20_GTonly.vcf -I

## but some GT=0/0 were missing, not sure why. sorted with grep, see rna_chr22.sh


## phaser: re-use phaser.array.sh with GBR22.sh

#####
## Prepare inputs for eQTL analysis
####################################

################################################################
## extract GT information from each sample and merge it with ASE information

vcf4AS /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/RNA 22 22

## format in inputs.R: save samples in "...for.AS.tab"
## convert tab files into vcf
tab2vcf /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/RNA for.AS.tab

## merge all samples into 1 vcf per chr

mergevcf /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE/RNA 22 22

## rm intermediate files: not done


## QC vcf file, make sure GT is in the right format: inputs.R
## make sure counts and vcf have the same samples:inputs.R

################## Run model for chr22 ##########

## Some analysis in comp.out.chr22.GT.noGT.lm.R and more updated in out.chr22.lm.GT.noGT.rna.v2.R


###########################################################################################
##### Effect of reference panel on phasing  ###

## select samples for : MXL

## good panel: EUR + AMR except MXL

## bad panel: AFR + South Asian

## for samples: take chr22 and make a vcf unphased.


## run in pop_RefPanel.R
