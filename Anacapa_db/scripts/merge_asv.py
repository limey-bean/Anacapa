# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program

# example usage python merge_asv.py nochim_forwardCO1.txt nochim_mergedCO1.txt nochim_reverseCO1.txt nochim_unmergedCO1.txt -o nochimCO1.txt

import argparse
import pandas as pd
import os

def read_table(file_name, rename=None):
    df = pd.read_csv(file_name, sep='\t', header=0, dtype='object')
    columns = list(df)
    old_id = columns[0]
    if rename:
        new_id = rename
    else:
        new_id = '_'.join(old_id.split('_')[-3:])

    df = df.rename(columns={old_id: new_id})
    sequence_column_names = ['sequence', 'sequencesF', 'sequencesR']
    for name in sequence_column_names:
        if name in columns:
            df = df.drop([name], axis=1)

    return df


def file_has_content(file_name):
    try:
        size = os.path.getsize(file_name)
        return size > 0
    except OSError:
        return False


def merge(file_names, outfile_name, rename=None):
    file_names = [fn for fn in file_names if file_has_content(fn)]
    main_df = read_table(file_names[0], rename)
    columns = list(main_df)
    for file_name in file_names[1:]:
        to_append = read_table(file_name, rename)
        main_df = main_df.append(to_append, ignore_index=True).fillna('0')
    all_columns = columns + [c for c in list(main_df) if c not in columns]
    main_df = main_df.loc[:, all_columns]
    main_df.to_csv(outfile_name, sep='\t', index=False)


parser = argparse.ArgumentParser(description='Merges ASV tables to a single one')
parser.add_argument('-o', '--out', type=str,
                    help='Write to a different output file (by default writes to merged.txt)')
parser.add_argument('-i', '--id', type=str,
                    help='rename the seq number column to specified name')
parser.add_argument('files', nargs='+',  help='Files to merge')

if __name__ == "__main__":
    args = parser.parse_args()
    files = args.files

    if args.out:
        outfile = args.out
    else:
        outfile = 'merged.txt'

    merge(files, outfile, rename=args.id)

