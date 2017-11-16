# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program
# python ~/group_alignments_to_files.py <outdir> <sam> <percent to keep> <allowable overhang> -p (if paired)


import re
import sys
import os
import argparse
import shutil
from collections import defaultdict, Counter


class SamEntry(object):
    def __init__(self, raw_row):
        self.qname = raw_row[0]
        self.flag = raw_row[1]
        self.rname = raw_row[2]
        self.pos = raw_row[3]
        self.mapq = raw_row[4]
        self.cigar = raw_row[5]
        self.rnext = raw_row[6]
        self.pnext = raw_row[7]
        self.tlen = raw_row[8]
        self.seq = raw_row[9]
        self.qual = raw_row[10]
        self.alignment_scores = [int(score.split(':')[-1]) for score in raw_row[11:13]]
        # sometimes it is a column early  because there is no second best match
        if raw_row[17][:2] == 'MD':
            self.md_z_flags = raw_row[17].split(':')[-1]
        else:
            self.md_z_flags = raw_row[18].split(':')[-1]

        total_match_count, total_count = self.calulate_match_count(self.md_z_flags)
        self.identity_ratio = total_match_count / total_count
        self.total_match = total_match_count

    def calulate_match_count(self, md_z_flags):
        tokenizer = re.compile(r'(\d+)|(\^[A-Z])|([A-Z])')
        total_match_count = 0.0
        total_mismatch_count = 0.0
        for item in tokenizer.finditer(md_z_flags):
            match_count, deletion, mismatch = item.groups()
            if match_count:
                total_match_count += int(match_count)
            if deletion:
                total_mismatch_count += 1
            if mismatch:
                total_mismatch_count += 1

        return total_match_count, (total_match_count + total_mismatch_count)

    # I am not completely sure if this is the logic that you want to use for checking unmapped at the ends
    # This function returns the max S at either end of the CIGAR score, 0  if there is no S at both ends
    def cigar_max_s(self):
        tokenizer = re.compile(r'(?:\d+)|(?:[A-Z=])')
        cigar_elements = tokenizer.findall(self.cigar)

        if cigar_elements[1] == 'S':
            beginning_s = int(cigar_elements[0])
        else:
            beginning_s = 0

        if cigar_elements[-1] == 'S':
            ending_s = int(cigar_elements[-2])
        else:
            ending_s = 0

        return max(beginning_s, ending_s)

class SummaryAppender(object):
    def __init__(self, local_sam_name, end_to_end_sam_name, taxonomy_file_name):
        with open(taxonomy_file_name) as taxonomy_file:
            self.taxonomy = self.load_taxonomy_file(taxonomy_file)
        self.parsed_bowtie_output = {}
        self.cutoffs = [(0.99, ".99-1.0"), (0.97, ".97-.9899"), (0.95, ".95-.9699"), (0.90, ".90-.9499"),
                        (0.85, ".85-.8999"), (0.80, ".80-.8499"), (0.00, ".00-.7999")]

        with open(local_sam_name) as sam_file:
            self.parsed_bowtie_output.update(self.parse_bowtie_output(sam_file, 'local'))
        with open(end_to_end_sam_name) as sam_file:
            self.parsed_bowtie_output.update(self.parse_bowtie_output(sam_file, 'end_to_end'))

    def load_taxonomy_file(self, taxonomy_file):
        taxonomy = {}
        for line in taxonomy_file:
            identifier = line.split('\t')[0].strip()
            path = line.split('\t')[1]
            taxonomy[identifier] = path.strip()
        return taxonomy

    def parse_bowtie_output(self, sam_file, mode):
        summary = {}
        previous_entries = []

        for line in sam_file:
            pieces = line.strip().split('\t')
            entry = SamEntry(pieces)
            if len(previous_entries) == 0 or previous_entries[0].qname == entry.qname:
                previous_entries.append(entry)
            else:
                summary[previous_entries[0].qname] = self.handle_group_of_entries(previous_entries, mode)
                previous_entries = [entry]

        # need to handle the stragglers
        if len(previous_entries) > 0:
            summary[previous_entries[0].qname] = self.handle_group_of_entries(previous_entries, mode)

        return summary

    def handle_group_of_entries(self, entries, mode):
        ids = [entry.identity_ratio for entry in entries]
        overhangs = [entry.cigar_max_s() for entry in entries]
        id_cutoff, id_group_name = self.calculate_cutoff_group(max(ids))
        taxonomies_in_range = [self.taxonomy[entry.rname] for entry in entries if entry.identity_ratio >= id_cutoff]
        accession_in_range = [entry.rname for entry in entries if entry.identity_ratio >= id_cutoff]
        taxonomy_consensus = self.determine_taxonomy_consensus(taxonomies_in_range)
        match_lengths = [entry.total_match for entry in entries]
        group_info = {}
        #, , taxonomic_paths_for_all_hits_within_percent_id_bin
        group_info['percent_ids'] = '{:0.4f}-{:0.4f}'.format(max(ids), min(ids))
        group_info['percent_id_bin'] = id_group_name
        group_info['single_or_multiple_hit'] = 'single' if len(entries) == 1 else 'multiple'
        group_info['length_of_hit'] = '{}-{}'.format(max(match_lengths), min(match_lengths))
        group_info['overhang_or_soft_clipping'] = '{}-{}'.format(max(overhangs), min(overhangs))
        group_info['consensus_taxonomic_path'] = taxonomy_consensus if taxonomy_consensus else 'not_assignable'
        group_info['taxonomic_paths_above_cutoff'] = '|'.join(taxonomies_in_range)
        group_info['accession_above_cutoff'] =  '|'.join(accession_in_range)
        group_info['mode'] = mode
        return group_info

    def calculate_cutoff_group(self, identity_ratio):
        for cutoff_value, possible_group_name in self.cutoffs:
            if identity_ratio >= cutoff_value:
                group_name = possible_group_name
                break
        return cutoff_value, group_name

    def determine_taxonomy_consensus(self, taxonomies):
        if len(taxonomies) == 1:
            return taxonomies[0]
        split_taxonomies = [t.split(';') for t in taxonomies]
        common_path = os.path.commonprefix(split_taxonomies)
        if len(common_path) == 0:
            return None

        return ';'.join(common_path)

    def append_to_summary(self, summary_file_name):
        fields = ['percent_ids','percent_id_bin','mode', 'single_or_multiple_hit', 'length_of_hit', 'overhang_or_soft_clipping',
                  'consensus_taxonomic_path', 'taxonomic_paths_above_cutoff', 'accession_above_cutoff']
        no_hit_entry = ['', '', '', '', '', '', 'no_hit', 'no_hit', '']

        current_summary_file = open(summary_file_name).readlines()
        current_header = current_summary_file[0]
        new_summary_file = open(summary_file_name + '.tmp', 'w')
        new_summary_file.write('\t'.join([current_header.rstrip('\n')] + fields) + '\n')
        for line in current_summary_file[1:]:
            name = line.split('\t')[0]
            old_line = line.rstrip('\n')
            if name in self.parsed_bowtie_output:
                sam_info = self.parsed_bowtie_output[name]
                ordered_sam_info = [sam_info[field] for field in fields]
                new_summary_file.write('\t'.join([old_line] + ordered_sam_info) + '\n')
            else:
                new_summary_file.write('\t'.join([old_line] + no_hit_entry))
        new_summary_file.close()
        shutil.move(summary_file_name + '.tmp', summary_file_name)

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


parser = argparse.ArgumentParser(description='Creates a summary table directly from bowtie2 output')
parser.add_argument('-v', '--verbose', help='Print out debugging information',
                    action='store_true')

parser.add_argument('summary_file', type=str,  help='File to append to')
parser.add_argument('taxonomy_file', type=str, help='Taxonomy file')
parser.add_argument('local_sam',  type=str, help='The local mode SAM file')
parser.add_argument('end_to_end_sam', type=str, help='The end-to-end SAM file')

if __name__ == '__main__':
    args = parser.parse_args()

    local_sam = args.end_to_end_sam
    end_to_end_sam = args.end_to_end_sam
    summary_file_name = args.summary_file
    taxonomy_file_name = args.taxonomy_file
    appender = SummaryAppender(local_sam, end_to_end_sam, taxonomy_file_name)
    appender.append_to_summary(summary_file_name)
