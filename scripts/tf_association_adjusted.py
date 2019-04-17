"""
File: tf_association_adjusted.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Run differential expression while adjusting for TF expression
"""
import pandas as pd
import argparse
import tf_association as tfa


def parse_arguments():
    info = 'Run differential expression while adjusting for TF expression'
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
    parser.add_argument('-l', '--log-transform',
                        action='store_true', default=False,
                        help='Whether to log transform the expression')
    #parser.add_argument('-c', '--copy-number',
                        #type=str, required=True,
                        #help='Copy number data')
    parser.add_argument('-a', '--adjusted-tfs',
                        type=str, required=True,
                        help='List of transcription factors for adjustment')
    parser.add_argument('-g', '--gene-name',
                        type=str, required=True,
                        help='Gene name of analysis')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read top TFs
    with open(opts['adjusted_tfs']) as handle:
        tfs = [l.strip() for l in handle]

    # output path
    output_path = opts['output']

    # analyze each TF for adjustment
    opts['output'] = None
    rename_dict = {
        'KMT2A': 'MLL',
        'KMT2B': 'MLL2',
        'KMT2C': 'MLL3',
        'KMT2D': 'MLL4',
        'NELFE': 'RDBP',
        'NELFA': 'WHSC2',
        'H2AZ': 'H2AFZ',
        'CHAMP1': 'ZNF828',
        'YBX3': 'CSDA',
        'HNRNPLL': 'HNRPLL',
        'STAT5': 'STAT5A',
        'ZBTB14': 'ZFP161',
        'ZFP69B': 'ZNF643'
    }
    result_list = []
    for tf in tfs:
        # some names are not their current hugo name
        if tf in rename_dict:
            tmp_tf = rename_dict[tf]
        else:
            tmp_tf = tf

        # run analysis
        result = tfa.main(opts,
                          gene_list=[opts['gene_name']],
                          covariate=tmp_tf)

        # format data
        result = result.rename(columns={opts['gene_name']: tf}).set_index('gene')
        result_list.append(result)

    # prepare output
    output_df = pd.concat(result_list, axis=1)
    output_df.to_csv(output_path, sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

