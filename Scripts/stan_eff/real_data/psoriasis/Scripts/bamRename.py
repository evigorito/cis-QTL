import os

def bam_rename(bamPath, outNames, outPath):
    """ Rule to created files with full path to bam and sample name and second file with full path to bam"""
    f=open(outNames, "w+")
    f2=open(outPath, "w+")
    for bam in bamPath:
        base=os.path.basename(os.path.dirname(bam))
        f.write(bam + " " + base + "\n")
        f2.write(bam + "\n")
    f.close()
    f2.close()

bam_rename(snakemake.input['bam'],
           str(snakemake.output['bamSampNames']),
           str(snakemake.output['bamPath']))


    
