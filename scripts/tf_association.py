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
import utils

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
    parser.add_argument('-im', '--immune',
                        type=str, required=False, default=None,
                        help='TCGA panimmune estimates')
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


def add_tumor_purity(expr_df, purity_path):
    """Add tumor purity to expression data"""
    # read in purity
    purity_df = pd.read_table(purity_path)
    purity_df['PatientID'] = purity_df['sample'].str[:12]
    purity_df = purity_df.drop_duplicates('PatientID')
    purity_df = purity_df.set_index('PatientID')

    # fill in purity
    expr_df['purity'] = purity_df['purity']
    expr_df['purity'] = expr_df['purity'].fillna(expr_df['purity'].mean())

    return expr_df


def add_tumor_subtype(expr_df, subtype_path, is_pancan=False):
    """Add subtype information to expression data."""
    # read in subtype information
    subtype_df = pd.read_table(subtype_path)

    # relate patient ID to subtype
    if is_pancan:
        subtype_df['SUBTYPE'] = subtype_df['DISEASE'] + '_' + subtype_df['SUBTYPE']
    patient2subtype = subtype_df.groupby('PATIENT_BARCODE')['SUBTYPE'].first().to_dict()

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
                       left_index=True, right_index=True,
                       how='left')

    return expr_df, unique_subtypes


def add_leukocyte_infiltration(expr_df, immune_path):
    # read in immune estimates
    immune_df = pd.read_table(immune_path)
    immune_df = immune_df.set_index('TCGA Participant Barcode')

    # merge in leukoctye fraction
    mycols = ['Leukocyte Fraction']
    expr_df = pd.merge(expr_df, immune_df[mycols],
                       left_index=True, right_index=True,
                       how='left')

    # fill missing values
    expr_df['Leukocyte Fraction'] = expr_df['Leukocyte Fraction'].fillna(expr_df['Leukocyte Fraction'].mean())

    return expr_df


def main(opts, gene_list=None, covariate=None):
    # read in mutations
    mut_df = pd.read_table(opts['mutation'])

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

    # add tumor purity
    if opts['tumor_purity']:
        expr_df = add_tumor_purity(expr_df, opts['tumor_purity'])

    # add immune infiltration
    if opts['immune']:
        expr_df = add_leukocyte_infiltration(expr_df, opts['immune'])

    # add subtype information
    is_pancan = os.path.basename(opts['mutation'])[:-4] == 'PANCAN'
    expr_df, unique_subtypes = add_tumor_subtype(expr_df, opts['subtype'], is_pancan=is_pancan)
    num_subtypes = len(unique_subtypes)

    # iterate through each gene
    output_dict = {}
    for gene in gene_list:
        # if pancan, check ctypes
        if is_pancan:
            mut_dir = os.path.dirname(opts['mutation'])
            ctypes = utils.check_signif_cancer_types(mut_dir, gene)
            if len(ctypes) == 1: continue
            is_in_ctypes = expr_df['subtype'].str.split('_', expand=True)[0].str[8:].isin(ctypes)
            expr_df2 = expr_df[is_in_ctypes].copy()
        else:
            expr_df2 = expr_df

        # identify mut vs wt samps
        mut_samps = mut_df.loc[mut_df[gene]==1, 'PatientID'].tolist()
        wt_samps = mut_df.loc[mut_df[gene]==0, 'PatientID'].tolist()

        # label mutations
        expr_df2['mutation'] = 0
        is_mut_samp = expr_df2.index.isin(mut_samps)
        expr_df2.loc[is_mut_samp, 'mutation'] = 1

        # create feature matrix
        if num_subtypes > 1:
            # remove any column which doesn't have at least two unique values
            unique_subtypes_cpy  = unique_subtypes.copy()
            num_uniq = expr_df2.nunique()
            constant_columns = num_uniq[num_uniq<2].index.tolist()
            for c in (set(constant_columns) & set(unique_subtypes_cpy)):
                unique_subtypes_cpy.remove(c)

            # now only add flags for the relevant subtypes
            mycols = unique_subtypes_cpy[:-1] + ['intercept', 'mutation']
        else:
            mycols = ['intercept', 'mutation']
        if opts['tumor_purity']:
            mycols += ['purity']
        if opts['immune'] and expr_df2['Leukocyte Fraction'].isnull().sum()==0:
            mycols += ['Leukocyte Fraction']

        # add additional covariate if specified
        if covariate:
            mycols.append(covariate)
        X = expr_df2[mycols].copy()
        if covariate:
            X[covariate] = X[covariate].fillna(X[covariate].mean())

        # drop any missing
        X = X.dropna()

        if len(X):
            # regress against all genes
            output_list = []
            for gene2 in (set(all_genes) - set([gene])):
                y = expr_df2[gene2]
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
