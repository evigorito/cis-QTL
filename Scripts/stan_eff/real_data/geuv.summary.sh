########################################################
## Apply bayesian trecase to GEUV data of EUR ancestry
########################################################

module load samtools/1.4.1

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


################################### Extract snp coordinates with REF and ALT alleles

vcf4rasqual /mrc-bsu/scratch/ev250/EGEUV1/quant/ASE 22 22 ASE.allsamples.vcf.gz  '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' /mrc-bsu/scratch/ev250/EGEUV1/quant/Btrecase/input snp.coord.txt

################# Prepare counts per gene and format inputs for Btrecase #######

## done in /home/ev250/Bayesian_inf/trecase/Scripts/stan_eff/real_data/QTL_input/inputs.R

################## Run model for chr22 ##########

## Started in 'top.snp.set.22.R' for gene-snp associations based on Chris eqtl geuvadis data.
## In this script I call others to run the top snp but also a snp in low LD to set rules to avoid running unecessary tests.

## Expanding into chr22.R, in inputs.R I could make a list of genes to test based on chr and level of expression so I can then make an array variable that goes through that list to run eQTL.

chr22.R
Error in dyn.load(libLFile) : 
  unable to load shared object '/tmp/RtmpvDL73J/file4ff67a20d.so':
  `maximal number of DLLs reached...

## options: remove libraries before running stan
## remove dll stan files as I go along: https://github.com/stan-dev/rstan/issues/448

dso_filename = mod@dso@dso_filename
  loaded_dlls = getLoadedDLLs()

dso_filename = model@dso@dso_filename
  loaded_dlls = getLoadedDLLs()
  if (dso_filename %in% names(loaded_dlls)) {
    message("Unloading DLL for model dso ", dso_filename)
    model.dll = loaded_dlls[[dso_filename]][['path']]
    dyn.unload(model.dll)
  } else {
    message("No loaded DLL for model dso ", dso_filename)
  }

  loaded_dlls = getLoadedDLLs()
  loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
  loaded_dlls <- loaded_dlls[grep("^file", names(loaded_dlls),value=T)]
  if (length(loaded_dlls) > 10) {
    for (dll in head(loaded_dlls, -10)) {
      message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
      dyn.unload(dll[['path']])
    }
  }
  message("DLL Count = ", length(getLoadedDLLs()), ": [", str_c(names(loaded_dlls), collapse = ","), "]")
