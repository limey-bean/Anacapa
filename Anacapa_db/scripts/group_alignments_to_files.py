# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program
# python ~/group_alignments_to_files.py <outdir> <sam> <percent to keep> <allowable overhang> -p (if paired)


import re
import sys
import argparse


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
    def should_keep_or_reject(self, entry):
        '''This function decides where to keep or reject a sequence'''

        if entry.cigar_max_s() > self.max_allowable_cigar_s:
            return False, 'not_mapped_at_ends'

        elif entry.alignment_scores[0] <= entry.alignment_scores[1]:
            return False, 'multiple_hits'

        elif entry.identity_ratio < self.identity_cutoff_to_keep:
            return False, 'low_percent_id'

        else:
            return True, None

    def choose_to_keep_or_reject(self, entry):
        should_keep, reason = self.should_keep_or_reject(entry)
        if should_keep:
            if self.verbose:
                print 'Writing {} to good file'.format(entry.qname)
            self.send_to_good_file(entry)
        else:
            if self.verbose:
                print 'Writing {} to reject file with reason {}'.format(entry.qname, reason)
            self.send_to_reject_file(reason, entry)

    def choose_to_keep_or_reject_pair(self, forward, backward):
        should_keep_forward, forward_reason = self.should_keep_or_reject(forward)
        should_keep_backward, backward_reason = self.should_keep_or_reject(backward)

        if should_keep_forward and should_keep_backward:
            if self.verbose:
                print 'Writing {}/{} to good file'.format(forward.qname, backward.qname)
            self.send_to_good_file(forward)
        if not should_keep_forward:
            if self.verbose:
                print 'Writing {}/{} to bad file with reason {}'.format(forward.qname, backward.qname,
                                                                        forward_reason)
            self.send_pair_to_reject_file(forward_reason, forward, backward)
        else:  # backward is the reason
            if self.verbose:
                print 'Writing {}/{} to bad file with reason {}'.format(forward.qname, backward.qname,
                                                                        backward_reason)
            self.send_pair_to_reject_file(backward_reason, forward, backward)

    def __init__(self, directory, max_allowable_cigar_s, identity_cutoff_to_keep, good_file_name='bowtie2_good_hits.txt',
                 general_reject_prefix='bowtie2_rejects_', verbose=False):

        self.directory = directory
        self.good_file_name = good_file_name
        self.cutoffs = [(0.99, ".99-1.0"), (0.97, ".97-.9899"), (0.95, ".95-.9699"), (0.90, ".90-.9499"),
                        (0.85, ".85-.8999"), (0.80, ".80-.8499"), (0.00, ".00-.7999")]
        self.general_reject_prefix = general_reject_prefix
        self.max_allowable_cigar_s = max_allowable_cigar_s
        self.identity_cutoff_to_keep = float(identity_cutoff_to_keep)
        self.verbose = verbose
        # keep all files open until the end because opening and closing is very slow
        self.file_cache = {}

    def send_pair_to_reject_file(self, prefix, forward, reverse):
        min_ratio = min(forward.identity_ratio, reverse.identity_ratio)
        group_name = self.calculate_cutoff_group(min_ratio)
        self.send_to_reject_file(prefix, forward, postfix='_forward', override_group_name=group_name)
        self.send_to_reject_file(prefix, reverse, postfix='_reverse', override_group_name=group_name)

    def calculate_cutoff_group(self, identity_ratio):
        for cutoff_value, possible_group_name in self.cutoffs:
            if identity_ratio >= cutoff_value:
                group_name = possible_group_name
                break
        return group_name

    def send_to_reject_file(self, prefix, entry, postfix='', override_group_name=None):
        # calculate cutoff

        if override_group_name:
            group_name = override_group_name
        else:
            group_name = self.calculate_cutoff_group(entry.identity_ratio)


        file_name = "{}{}{}{}{}.fasta".format(self.directory, self.general_reject_prefix,
                                              prefix, group_name, postfix)

        content = ">{}\n{}\n".format(entry.qname, entry.seq)

        self.write_with_cache(file_name, content)

    def send_to_good_file(self, entry):
        file_name = self.directory + self.good_file_name
        content = "{}\t{}\t{:0.4f}\n".format(entry.qname, entry.rname, entry.identity_ratio)
        self.write_with_cache(file_name, content)

    def clean_up_file_cache(self):
        for name, f in self.file_cache.items():
            if self.verbose:
                print 'Saving file {}'.format(name)
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

    def process_paired_sam_file(self, sam_file_name):
        with open(sam_file_name) as sam_file:
            while True:
                line1 = sam_file.readline().strip()
                line2 = sam_file.readline().strip()
                if not line2:
                    break
                forward = SamEntry(line1.split('\t'))
                backward = SamEntry(line2.split('\t'))
                self.choose_to_keep_or_reject_pair(forward, backward)
        self.clean_up_file_cache()


def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


parser = argparse.ArgumentParser(description='Groups output of bowtie2 into files')
parser.add_argument('-v', '--verbose', help='Print out debugging information',
                    action='store_true')
parser.add_argument('-p', '--paired', help='Paired mode, input file is treated as '
                                           'pairs of forward and backward sequences',
                    action='store_true')
parser.add_argument('output_directory', type=str, help='Directory where output will be written')
parser.add_argument('input_file', type=str, help='The SAM file to be read')

parser.add_argument('s_to_allow', type=float, help='Allowable query overhang')

parser.add_argument('keep_percent', type=restricted_float, help='percent to keep')

if __name__ == '__main__':
    args = parser.parse_args()

    output_directory = args.output_directory
    if not output_directory.endswith('/'):
        output_directory = output_directory + '/'

    sam_file_name = args.input_file
    overhang = args.s_to_allow
    keep = args.keep_percent

    bowtie_sorter = BowtieSorter(output_directory, overhang, keep, verbose=args.verbose)
    if args.paired:
        bowtie_sorter.process_paired_sam_file(sam_file_name)
    else:
        bowtie_sorter.process_sam_file(sam_file_name)

