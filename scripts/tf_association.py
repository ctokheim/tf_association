"""
File: diff_expr.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Analyze differential expression
"""
import argparse
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import statsmodels.api as sm
from numpy.random import RandomState
import csv
import os

import warnings
warnings.filterwarnings("error")


def parse_arguments():
    info = 'Analyze differential expression'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Expression data')
    parser.add_argument('-s', '--subtype',
                        type=str, required=True,
                        help='Cancer subtype information')
    parser.add_argument('-t', '--tumor-purity',
                        type=str, required=False, default=None,
                        help='Tumor purity')
    parser.add_argument('-m', '--mutation',
                        type=str, required=True,
                        help='Mutation data')
    #parser.add_argument('-c', '--copy-number',
                        #type=str, required=True,
                        #help='Copy number data')
    parser.add_argument('-l', '--log-transform',
                        action='store_true', default=False,
                        help='Whether to log transform the expression')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts, gene_list=None, covariate=None):
    # read in mutations
    mut_df = pd.read_table(opts['mutation'])

    # read in subtype information
    subtype_df = pd.read_table(opts['subtype'])
    patient2subtype = subtype_df.groupby('PATIENT_BARCODE')['SUBTYPE'].first().to_dict()

    # read expression data
    expr_df = pd.read_table(opts['input'])
    expr_df = expr_df.drop_duplicates(subset=['gene_id'])
    expr_df = expr_df.set_index('gene_id').T
    if opts['log_transform']:
        expr_df = np.log2(expr_df+1)

    # complete list of genes with available expression data
    if not gene_list:
        gene_list = mut_df.columns[1:].tolist()
    all_genes = list(set(expr_df.columns.tolist()))

    # add intercept term
    expr_df['intercept'] = 1

    if opts['tumor_purity']:
        # read in purity
        purity_df = pd.read_table(opts['tumor_purity'])
        purity_df['PatientID'] = purity_df['sample'].str[:12]
        purity_df = purity_df.drop_duplicates('PatientID')
        purity_df = purity_df.set_index('PatientID')

        # fill in purity
        expr_df['purity'] = purity_df['purity']
        expr_df['purity'] = expr_df['purity'].fillna(expr_df['purity'].mean())

    # fill in tumor subtypes
    expr_df['subtype'] = 'Subtype:' + expr_df.index.map(patient2subtype).fillna('Not_Applicable')
    unique_subtypes = list(expr_df['subtype'].unique())
    num_subtypes = len(unique_subtypes)
    tmp = expr_df.reset_index().rename(columns={'index':'PatientID'})[['intercept', 'PatientID', 'subtype']]
    subtype_onehot = (
        tmp.pivot(values='intercept',
                  index='PatientID',
                  columns='subtype')
                  .fillna(0)
                  .astype(int)
    )
    expr_df = pd.merge(expr_df, subtype_onehot,
                       left_index=True, right_index=True, how='left')

    # iterate through each gene
    output_dict = {}
    for gene in gene_list:
        mut_samps = mut_df.loc[mut_df[gene]==1, 'PatientID'].tolist()
        wt_samps = mut_df.loc[mut_df[gene]==0, 'PatientID'].tolist()

        # label mutations
        expr_df['mutation'] = 0
        is_mut_samp = expr_df.index.isin(mut_samps)
        expr_df.loc[is_mut_samp, 'mutation'] = 1

        # create feature matrix
        if num_subtypes > 1:
            mycols = unique_subtypes[:-1] + ['intercept', 'mutation']
        else:
            mycols = ['intercept', 'mutation']
        if opts['tumor_purity']:
            mycols += ['purity']
        # add additional covariate if specified
        if covariate:
            mycols.append(covariate)
        X = expr_df[mycols].copy()
        if covariate:
            X[covariate] = X[covariate].fillna(X[covariate].mean())

        # regress against all genes
        output_list = []
        for gene2 in (set(all_genes) - set([gene])):
            y = expr_df[gene2]
            result = sm.OLS(y, X).fit()
            t_stats = result.params / result.bse
            mut_t_stat = t_stats.loc['mutation']
            output_list.append([gene2, mut_t_stat])

        # save output
        output_series = pd.DataFrame(output_list, columns=['gene', gene]).set_index('gene')[gene]
        output_dict[gene] = output_series

    # save output
    rename_dict = {'index': 'gene'}
    output_df = pd.DataFrame(output_dict).reset_index().rename(columns=rename_dict)
    if opts['output']: output_df.to_csv(opts['output'], sep='\t', index=False)

    return output_df


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
