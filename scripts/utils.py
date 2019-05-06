import csv
import pandas as pd
import os
import glob
import csv


def read_rabit(path):
    """Read the output file created by Rabit."""
    # read data
    with open(path) as handle:
        myreader = csv.reader(handle, delimiter='\t')

        # iterate through each line
        result_list = []
        list_of_lines = []
        for line in myreader:
            # if header then add to result list
            is_header = line[0].startswith('>')
            if is_header and list_of_lines:
                result_list.append(list_of_lines)
                list_of_lines= []

            # add line
            list_of_lines.append(line)
        result_list.append(list_of_lines)

    # return None if empty
    if result_list == [[]]:
        return None

    # convert to dataframe
    df_list = [read_rabit_section(l) for l in result_list]

    return df_list


def read_rabit_section(mylist):
    """Read a section within the rabit output file."""
    # create dataframe
    df = pd.DataFrame(mylist[1:], columns=mylist[0])

    # convert dtype
    for c in df.columns[1:]:
        df[c] = df[c].astype(float)

    # add gene name
    gene = mylist[0][0][1:]
    df['gene'] = gene
    # rename first column to avoid idiosyncratic behavior
    df = df.rename(columns={'>'+gene: 'chipseq'})

    # drop non chip-seq rows
    is_non_chipseq = ~df['chipseq'].isin(['Degree', 'CpG', 'Intercept'])
    df = df[is_non_chipseq]

    # add a column for the chipseq TF name
    df['tf'] = df.chipseq.str.extract('[0-9]+\.(.+)@[0-9]+')

    return df


def check_signif_cancer_types(mut_dir, gene):
    """Check which cancer types the gene is significant in"""
    # figure out which cancer types it is significant in
    ctype_list = []
    for f in glob.glob(os.path.join(mut_dir, '*.txt')):
        with open(f) as handle:
            myreader = csv.reader(handle, delimiter='\t')
            header = next(myreader)
            if gene in header:
                ctype = os.path.basename(f)[:-4]
                if ctype != 'PANCAN':
                    ctype_list.append(ctype)

    return ctype_list
