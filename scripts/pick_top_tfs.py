"""
File: pick_top_tfs.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Pick out the most significant TFs
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
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output directory')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df_list = []
    for f in glob.glob(os.path.join(opts['input'], '*.txt')):
        # read data
        tmp_df_list = utils.read_rabit(f)

        if tmp_df_list:
            tmp_df = pd.concat(tmp_df_list)

            # add cancer type
            ctype = os.path.basename(f)[:-4]
            tmp_df['Cancer Type'] = ctype
            df_list.append(tmp_df)

    # concatenate results
    df = pd.concat(df_list)
    df = df.sort_values('Pr(>|t|)')

    # pick top 10
    top10_df = df.groupby(['Cancer Type', 'gene']).head(10)
    top10_df = top10_df.sort_values(['Cancer Type', 'gene'])
    top10_df.index = range(len(top10_df))

    # extract gene name of top tf
    top10_df['tf'] = top10_df.chipseq.str.extract('[0-9]+\.(.+)@[0-9]+')

    # iterate through each group
    index_dict = top10_df.groupby(['gene', 'Cancer Type']).indices
    for gene, ctype in index_dict:
        ixs = index_dict[(gene, ctype)]
        outpath = os.path.join(opts['output'], '{0}.{1}.txt'.format(gene, ctype))
        top10_df.loc[ixs, 'tf'].to_csv(outpath, sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


