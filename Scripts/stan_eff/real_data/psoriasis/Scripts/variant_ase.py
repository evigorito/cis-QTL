import pandas;
import math;
import sys;

def main(haplotypic_counts, out, id_separator="_",  min_cov=0):
    """ This function assigns ASE counts to variants avoiding double counting  when one read overlaps more than 1 variant. 
    Inputs: 
    ## required
    'haplotypic_counts' file, output file from phASER containing read counts for haplotype blocks. NOTE: unphased_vars must have been enabled when phASER was run.
    'out', full name for output file
    ## optional, phaser options to make it compatible with phaser run
    "id_separator", default="_", help="Separator used for generating unique variant IDs when phASER was run."
    "min_cov", type=int, default=0, help="Minimum total coverage for a variant to be outputted."
    
    Output:
    File with: contig position variantID refAllele altAllele refCounts, altCounts and totalCounts for each het variant. Each file corresponds to one bam file. 

"""

    df_haplo_counts_master = pandas.read_csv(haplotypic_counts, sep="\t", index_col=False);

    if "bam" not in df_haplo_counts_master.columns:
        print("ERROR - this function is only compatible with results from phASER v1.0.0+");
        sys.exit(1)

    stream_out = open(out, "w");
    stream_out.write("\t".join(["contig","position","variantID","refAllele","altAllele","refCount","altCount","totalCount","bam"])+"\n");

    # keep record of variants
    dict_variants = {};

    # produce a separate output file for each bam
    for xbam in set(df_haplo_counts_master['bam']):
        df_haplo_counts = df_haplo_counts_master[(df_haplo_counts_master.bam == xbam)];
        for index, row in df_haplo_counts.iterrows():
            chrom = str(row['contig']);
            if row['totalCount'] > 0:
                xvars = row['variants'].split(",");
                if id_separator not in xvars[0] or xvars[0].count(id_separator) < 3:
                    print("ERROR - ID separator not found in variant ID, please ensure that id_separator is set correctly.")
                    sys.exit(1)
                    
                # collect used reads whithin the haplotype block to avoid double-counting
                hap_a_reads = []; 
                hap_b_reads = [];
                ##print('xvars ' + row['variants']  +  ' hap_a_reads ' + ",".join(str(e) for e in hap_a_reads))
                for xvar in xvars:
                    xvar_index = xvars.index(xvar);
                    fields = xvar.split(id_separator);
                    xvar_pos = int(fields[1]);
                    refAllele = fields[2];
                    altAllele  = fields[3];
                    hapA_allele = str(row['haplotypeA']).split(",")[xvar_index]

                    dict_variants[xvar] = {'contig':chrom,'position':xvar_pos,'refAllele':refAllele, 'altAllele':altAllele, 'refCount':0, 'altCount':0, 'totalCount':0, 'bam': xbam};
                    if len(xvars) == 1:
                        # if haplotype only has one variant then there is no double counting
                        aCount = row['aCount']
                        bCount = row['bCount']

                    else:
                        hap_a = str(row['aReads']).split(";")[xvar_index].split(",");
                        hap_b = str(row['bReads']).split(";")[xvar_index].split(",");
                        # remove blank read IDs (these are created when there are no reads mapping to a variant on a given haplotype
                        if "" in hap_a: hap_a.remove("");
                        if "" in hap_b: hap_b.remove("");                     
                        aCount = len(set(x for x in hap_a if x not in set(hap_a_reads)));
                        bCount = len(set(x for x in hap_b if x not in set(hap_b_reads)));
                        hap_a_reads += hap_a;
                        hap_b_reads += hap_b;

                    if aCount + bCount >= min_cov:
                        ## assign counts to ref and alternative alleles, invert if need be
                        if hapA_allele == refAllele:    
                            refCount = aCount
                            altCount = bCount
                        else:
                            refCount = bCount
                            altCount = aCount

                        dict_variants[xvar] = {'contig':chrom,'position':xvar_pos,'variantID': xvar, 'refAllele':refAllele, 'altAllele':altAllele,'refCount':refCount, 'altCount':altCount, 'totalCount': refCount + altCount, 'bam': xbam}

                        ## save output
                        stream_out.write("\t".join(list(map(str,dict_variants[xvar].values()))) + "\n")

    stream_out.close();


main(snakemake.input['hapCounts'], snakemake.output[0], snakemake.params['sep'], int(snakemake.params['min_cov']) )


#main("/mrc-bsu/scratch/ev250/psoriasis/phaser/chr1.SRR1146100.haplotypic_counts.txt", "/mrc-bsu/scratch/ev250/psoriasis/ASE/chr1.SRR1146100.corrected.haplotypic_counts.txt")
