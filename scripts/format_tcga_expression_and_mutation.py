"""
File: format_tcga_gene_expression.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Format TCGA gene expression matrix to be more amenable
"""
import csv
import argparse
import pandas as pd
import os


def parse_arguments():
    info = 'Format gene expression matrix to be more amenable'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='TCGA RNAseq gene expression matrix')
    parser.add_argument('-c', '--cancer-type',
                        type=str, required=True,
                        help='Cancer type')
    parser.add_argument('-m', '--mutation',
                        type=str, required=True,
                        help='Mutation data')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output directory')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in mutations
    mut_df = pd.read_table(opts['mutation'])
    ctype = opts['cancer_type']
    if ctype != 'PANCAN':
        mut_df = mut_df[mut_df['CODE']==ctype]
    # save mutation data
    mut_path = os.path.join(opts['output'], 'mutation/{0}.txt'.format(ctype))
    mut_df.to_csv(mut_path, sep='\t', index=False)
    # create list of patient IDs
    samp_list = mut_df['Tumor_Sample_Barcode'].str[:12].tolist()

    # iterate over each line in the file
    with open(opts['input']) as handle:
        myreader = csv.reader(handle, delimiter='\t')

        # parse header
        header = next(myreader)
        # parse the relevant cols
        header_ixs = []
        header_cols = []
        for i in range(len(header)):
            if i==0:
                header_ixs.append(i)
                header_cols.append(header[i])
                continue

            # add relevant samps
            samp = header[i][:12]
            if samp in samp_list:
                header_ixs.append(i)
                header_cols.append(samp)

        output = [header_cols]
        for line in myreader:
            gene_name = line[0].split('|')[0]
            if gene_name == '?':
                continue

            output.append([gene_name]+[line[ix] for ix in header_ixs[1:]])

    # write output
    expr_path = os.path.join(opts['output'], 'expr/{0}.txt'.format(ctype))
    with open(expr_path, 'w') as write_handle:
        mywriter = csv.writer(write_handle, delimiter='\t', lineterminator='\n')
        mywriter.writerows(output)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

