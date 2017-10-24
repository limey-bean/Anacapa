# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program

from collections import defaultdict, Counter
import argparse
import os

def load_taxonomy_file(taxonomy_file):
    taxonomy = {}
    for line in taxonomy_file:
        identifier = line.split('\t')[0].strip()
        species = line.split('\t')[1]
        taxonomy[identifier] = species.strip()
    return taxonomy


def count_species(hits_file, taxonomy):
    # This is of the form {sample_id: {species_name: count}}
    # for example {'sample1': {'Phascolarctos cinereus': 100}}
    species_counts = defaultdict(lambda: Counter())
    for line in hits_file:
        sequence_id, taxonomy_id, percent = line.strip().split('\t')
        sample_name = sequence_id.split('_:')[0]
        species_name = taxonomy[taxonomy_id]
        species_counts[sample_name][species_name] += 1

    return species_counts


def create_output(output_file, species_counts):
    all_samples = species_counts.keys()

    all_species = set()
    for sample in species_counts.values():
        all_species.update(sample.keys())

    # write the header
    output_file.write('#OTU ID\t{}\n'.format('\t'.join(all_samples)))

    for species_name in all_species:
        counts = []
        for sample_name in all_samples:
            counts.append(species_counts[sample_name][species_name])
        count_columns = '\t'.join(str(c) for c in counts)
        output_file.write('{}\t{}\n'.format(species_name, count_columns))


def create_summary(hits_file_name, taxonomy_file_name, output_file_name):
    with open(taxonomy_file_name) as taxonomy_file:
        taxonomy = load_taxonomy_file(taxonomy_file)

    with open(hits_file_name) as hits_file:
        species_counts = count_species(hits_file, taxonomy)

    with open(output_file_name, 'w') as output_file:
        create_output(output_file, species_counts)


parser = argparse.ArgumentParser(description='Summarizes good bowtie2 hits into species counts per sample')
parser.add_argument('taxonomy_file', type=str, help='Taxonomy file used to look up species name from an identifier')
parser.add_argument('good_hits_file', type=str, help='File to summarize')
parser.add_argument('output_file', type=str, help='File where output should be written')

if __name__ == '__main__':
    args = parser.parse_args()
    hits_file_name = args.good_hits_file
    taxonomy_file_name = args.taxonomy_file
    output_file_name = args.output_file

    create_summary(hits_file_name, taxonomy_file_name, output_file_name)