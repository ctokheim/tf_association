"""
File: format_rabit_results.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Format the output from Rabit runs
"""
import pandas as pd
import argparse
import utils
import glob
import os


def parse_arguments():
    info = 'Pick out the most significant TFs'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Directory of rabit output file')
    parser.add_argument('-b', '--baseline',
                        type=str, required=True,
                        help='Directory of rabit output files without TF adjustments')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    output_list = []
    for f in glob.glob(os.path.join(opts['input'], '*.txt')):
        # figure out gene name / cancer type
        gene, ctype, _ = os.path.basename(f).split('.')

        # read data
        tmp_df_list = utils.read_rabit(f)

        if tmp_df_list:
            for df in tmp_df_list:
                tf = df['gene'].iloc[0]
                is_adjusted_gene = df['gene']==df['tf']
                pval_series = df.loc[is_adjusted_gene, 'Pr(>|t|)']

                # get p-value if listed by rabit
                if len(pval_series):
                    pval = pval_series.iloc[0]
                else:
                    pval = 0.5
                output_list.append([gene, ctype, tf, pval])

    # compile results into a dataframe
    mycols = ['gene', 'cancer type', 'tf', 'pvalue']
    output_df = pd.DataFrame(output_list, columns=mycols)

    # grab top 10 TFs from original results
    output_list = []
    for f in glob.glob(os.path.join(opts['baseline'], '*.txt')):
        # read data
        tmp_df_list = utils.read_rabit(f)

        if tmp_df_list:
            ctype = os.path.basename(f)[:-4]
            for df in tmp_df_list:
                tmp = df.sort_values('Pr(>|t|)').head(10).copy()
                tmp['cancer type'] = ctype
                output_list.append(tmp[['gene', 'cancer type', 'tf', 'Pr(>|t|)']])

    # concatenate original rabit results
    original_df = pd.concat(output_list)
    original_df = original_df.rename(columns={'Pr(>|t|)': 'original pvalue'})

    # merge in original results
    output_df = pd.merge(output_df, original_df,
                         on=['gene', 'cancer type', 'tf'],
                         how='left')

    # save output
    output_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


