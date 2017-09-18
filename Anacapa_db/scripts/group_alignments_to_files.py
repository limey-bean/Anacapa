# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program


import re
import sys

import re
import sys

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
        self.md_z_flags = raw_row[18].split(':')[-1]
        self.identity_ratio = self.calculate_identity_ratio(self.md_z_flags)

    def calculate_identity_ratio(self, md_z_flags):
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

        return total_match_count / (total_match_count + total_mismatch_count)

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


class BowtieSorter(object):
    ##########################################
    # Modify this function to change criteria !
    def choose_to_keep_or_reject(self, entry):
        '''This function decides where to keep or reject a sequence'''

        if entry.cigar_max_s() > self.max_allowable_cigar_s:
            self.send_to_reject_file('not_mapped_at_ends', entry)

        elif entry.alignment_scores[0] <= entry.alignment_scores[1]:
            self.send_to_reject_file('multiple_hits', entry)

        elif entry.identity_ratio < self.identity_cutoff_to_keep:
            self.send_to_reject_file('low_percent_id', entry)

        else:
            self.send_to_good_file(entry)


    def __init__(self, directory, good_file_name='bowtie2_good_hits.txt',
                 general_reject_prefix='bowtie2_rejects_', max_allowable_cigar_s=25):

        self.directory = directory
        self.good_file_name = good_file_name
        self.cutoffs = [(0.97, ".97-1.0"), (0.95, ".95-.9699"), (0.9, ".90-.9499"),
                        (0.80, ".80-.8999"), (0.0, ".00-.7999")]
        self.general_reject_prefix = general_reject_prefix
        self.max_allowable_cigar_s = max_allowable_cigar_s
        self.identity_cutoff_to_keep = 0.97
        # keep all files open until the end because opening and closing is very slow
        self.file_cache = {}



    def send_to_reject_file(self, prefix, entry):
        # calculate cutoff
        for cutoff_value, possible_group_name in self.cutoffs:
            if entry.identity_ratio >= cutoff_value:
                group_name = possible_group_name
                break

        file_name = "{}{}{}{}.fasta".format(self.directory, self.general_reject_prefix,
                                      prefix, group_name)

        content = ">{}\n{}\n".format(entry.qname, entry.seq)

        self.write_with_cache(file_name, content)


    def send_to_good_file(self, entry):
        file_name = self.directory + self.good_file_name
        content = "{}\t{}\t{:0.4f}\n".format(entry.qname, entry.rname, entry.identity_ratio)
        self.write_with_cache(file_name, content)

    def clean_up_file_cache(self):
        for f in self.file_cache.values():
            f.close()
        self.file_cache = {}

    def write_with_cache(self, file_name, content):
        if file_name in self.file_cache:
            file_to_write = self.file_cache[file_name]
        else:
            file_to_write = open(file_name, 'a')
            self.file_cache[file_name] = file_to_write

        file_to_write.write(content)

    def process_sam_file(self, sam_file_name):
        with open(sam_file_name) as sam_file:
            for line in sam_file:
                pieces = line.strip().split('\t')
                entry = SamEntry(pieces)
                self.choose_to_keep_or_reject(entry)
        self.clean_up_file_cache()



if __name__ == 'main':
    # TODO: use argparse
    output_directory = sys.argv[1]
    sam_file_name = sys.argv[2]
    bowtie_sorter = BowtieSorter(output_directory)
    bowtie_sorter.process_sam_file(sam_file_name)






