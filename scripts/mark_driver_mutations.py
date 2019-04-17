"""
File: mark_driver_mutations.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Flag driver mutations
"""
import pandas as pd
import argparse


def parse_arguments():
    info = 'Flag driver mutations'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='MAF file')
    parser.add_argument('-d', '--driver',
                        type=str, required=True,
                        help='File containing driver genes')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    driver_df = pd.read_table(opts['driver'])
    mut_df = pd.read_table(opts['input'])
    mut_df['PatientID'] = mut_df['Tumor_Sample_Barcode'].str[:12]

    # find out the relevant type
    is_signif = driver_df['qvalue']<=0.05
    is_high_score = driver_df['score']>0.5
    driver_df = driver_df[is_signif & is_high_score]
    driver_df['info2'] = driver_df['info'].str[5:]
    gene_dict = dict(driver_df.groupby('gene')['info2'].unique())

    # list of variant types
    og_var_list = ['Missense_Mutation']
    tsg_var_list = [
        'Nonsense_Mutation', 'Frame_Shift_Del', 'Splice_Site',
        'Frame_Shift_Del', 'Frame_Shift_Ins', 'Translation_Start_Site',
        'Nonstop_Mutation'
    ]
    driver_var_list = og_var_list + tsg_var_list

    # iterate through each gene
    output_dict = {}
    for gene in gene_dict:
        driver_types = list(gene_dict[gene])

        # flag mutations
        is_gene = mut_df['Hugo_Symbol'] == gene
        mut_df[gene] = 0
        if 'oncogene' in driver_types:
            is_var = mut_df['Variant_Classification'].isin(og_var_list)
            mut_df.loc[is_gene & is_var, gene] = 1
        if 'tsg' in driver_types:
            is_var = mut_df['Variant_Classification'].isin(tsg_var_list)
            mut_df.loc[is_gene & is_var, gene] = 1
        if len(driver_types)==1 and driver_types[0]=='driver':
            is_var = mut_df['Variant_Classification'].isin(driver_var_list)
            mut_df.loc[is_gene & is_var, gene] = 1

        # flag by patient id
        patient_flags = mut_df.groupby('PatientID')[gene].max()

        # add to output
        output_dict[gene] = patient_flags

    # save results
    output_df = pd.DataFrame(output_dict)
    output_df.to_csv(opts['output'], sep='\t')


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)

