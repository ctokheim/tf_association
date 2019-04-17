"""
File: convert_to_entrez_id.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Convert hugo gene names to entrez IDs
"""
import pandas as pd
import numpy as np
import argparse


def parse_arguments():
    info = 'Convert hugo gene names to entrez IDs'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Differential expression file')
    parser.add_argument('-c', '--conversion',
                        type=str, required=True,
                        help='Gene name to Entrez ID file')
    parser.add_argument('-e', '--expression',
                        type=str, required=True,
                        help='Expression data')
    parser.add_argument('-eq', '--expressed-quantile',
                        type=float, required=True,
                        help='Expression quantile for threshold')
    parser.add_argument('-t', '--top-variability',
                        type=int, required=True,
                        help='Number of genes to select based on variability')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Differential expression file with entrez ids')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read input
    df = pd.read_table(opts['input'])

    # read in expression data and filter to relevant subset
    expr_df = pd.read_table(opts['expression'], index_col=0)
    expr_df = np.log2(expr_df+1)
    meds = expr_df.median(axis=1).rank(pct=True)
    top_expressed = meds[meds>opts['expressed_quantile']].index
    cv = expr_df.loc[top_expressed, :].std(axis=1) / expr_df.loc[top_expressed, :].mean(axis=1)
    cv = cv.sort_values(ascending=False)
    most_variable = cv.iloc[:opts['top_variability']].index.tolist()

    # keep only genes based on exprssion filters
    #df = df[df['gene'].isin(most_variable)]

    # read in gene conversion file
    gene2id = pd.read_table(opts['conversion'], dtype={'gene': str, 'entrez_id':str})
    rename_dict = dict(gene2id.values)

    # replace gene names
    df['gene'] = df['gene'].replace(rename_dict)

    # set gene name as index
    df = df.dropna() # newly added
    df = df.drop_duplicates(subset=['gene'])
    df = df.set_index('gene')

    # save output
    # previously dropna occurred just before CSV
    df.to_csv(opts['output'], sep='\t', index_label=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

