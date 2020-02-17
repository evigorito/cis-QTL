from tabula import read_pdf
import pandas as pd

def get_pdf(pdf, out, page=[1]):
    df = read_pdf(pdf, pages = page[0],multiple_tables=True)
    if len(page) > 1:
        for p in range(1,len(page)):
            df.append(read_pdf(pdf, page[p]), ignore_index=True)
    df.to_csv(out)

get_pdf(pdf = str(snakemake.input), page=snakemake.params['pages'], out=str(snakemake.output) )
