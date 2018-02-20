# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program

# command to run: python path/to/append_blca_to_summary.py -o optional/path/to/output/txtfile path/to/summary_dada2.txt path/to/blca.out

import argparse
import shutil
from collections import namedtuple

TaxonomyEntry = namedtuple('TaxonomyEntry', ['values', 'confidences', 'accessions'])


def read_taxonomy(taxonomy_file_name):
    taxonomy = {}
    with open(taxonomy_file_name) as taxonomy_file:
        for line in taxonomy_file:
            pieces = line.strip().split("\t")
            if len(pieces) == 2:
                taxonomy[pieces[0]] = TaxonomyEntry(pieces[1], '', '')
            else:
                taxonomy[pieces[0]] = TaxonomyEntry(pieces[1], pieces[2], pieces[3])
    return taxonomy


def append_to_summary(summary_file_name, taxonomy_file_name, output_file_name):
    taxonomy = read_taxonomy(taxonomy_file_name)
    fields = ['taxonomy', 'taxonomy_confidence', 'accessions']

    missing_taxonomy = TaxonomyEntry('not_found', '', '')

    current_summary_file = open(summary_file_name).readlines()
    current_header = current_summary_file[0]
    new_summary_file = open(output_file_name + '.tmp', 'w')
    new_summary_file.write('\t'.join([current_header.rstrip('\n')] + fields) + '\n')
    for line in current_summary_file[1:]:
        name = line.split('\t')[0]
        old_line = line.rstrip('\n')
        entry = taxonomy.get(name, missing_taxonomy)
        new_summary_file.write('\t'.join([old_line] + list(entry)) + '\n')
    new_summary_file.close()
    shutil.move(output_file_name + '.tmp', output_file_name)


parser = argparse.ArgumentParser(description='Creates a summary table from BLCA output')
parser.add_argument('-o', '--out', type=str,
                    help='Write to a different output file (by default writes directly to input)')
parser.add_argument('summary_file', type=str,  help='File to append to')
parser.add_argument('blca_taxonomy', type=str,  help='BLCA generated taxonomy file')

if __name__ == '__main__':
    args = parser.parse_args()
    summary_file = args.summary_file
    blca_taxonomy = args.blca_taxonomy
    if args.out:
        outfile = args.out
    else:
        outfile = summary_file
    append_to_summary(summary_file, blca_taxonomy, outfile)
