# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program
# python ~/group_alignments_to_files.py <outdir> <sam> <percent to keep> <allowable overhang> -p (if paired)


import re
import sys
import os
import glob
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
    def __init__(self, bowtie_out_directory):

        self.parsed_bowtie_output = {}
        for file_name in glob.glob(bowtie_out_directory + '*.sam'):
            with open(file_name) as sam_file:
                if '_local.sam' in file_name:
                    self.parsed_bowtie_output.update(self.parse_bowtie_output(sam_file, 'local'))
                elif '_end_to_end.sam' in file_name:
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
        group_info = {}
        # multiple_or_single_hit end_to_end_or_local input_sequence_length
        group_info['single_or_multiple_hit'] = 'single' if len(entries) == 1 else 'multiple'
        group_info['end_to_end_or_local'] = mode
        group_info['input_sequence_length'] = str(max([len(entry.seq) for entry in entries]))
        group_info['max_percent_id'] = str(max([entry.identity_ratio for entry in entries]))

        return group_info

    def append_to_summary(self, summary_file_name):
        fields = ['single_or_multiple_hit', 'end_to_end_or_local', 'max_percent_id', 'input_sequence_length']
        no_hit_entry = ['', '', '', 'no_hit']

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
                new_summary_file.write('\t'.join([old_line] + no_hit_entry) + '\n')
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
parser.add_argument('bowtie_out_directory',  type=str, help='The directory where all the SAM files are kept')

if __name__ == '__main__':
    args = parser.parse_args()
    bowtie_out_directory = args.bowtie_out_directory
    summary_file_name = args.summary_file
    appender = SummaryAppender(bowtie_out_directory)
    appender.append_to_summary(summary_file_name)
