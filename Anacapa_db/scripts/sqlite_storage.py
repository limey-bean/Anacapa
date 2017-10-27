# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program


import sqlalchemy as sql
from sqlalchemy import Column, Integer, Float, Text, Index
import sqlite3
import argparse

class BowtiePairStorage:
    def __init__(self, db_path):
        self.engine = sql.create_engine('sqlite:///' + db_path)
        self.conn = self.engine.connect()

        metadata = sql.MetaData()
        self.table = sql.Table('bowtie_pairs', metadata,
                               Column('forward_name', Text),
                               Column('reverse_name', Text),
                               Column('forward_sequence', Text),
                               Column('reverse_sequence', Text),
                               Column('identity_ratio', Float),
                               Column('max_s', Integer),
                               Column('single_best_match', Integer),
                               Column('reject_reason', Integer)
                               )
        metadata.create_all(self.engine, checkfirst=True)

    def store_pair(self, forward, reverse, reject_reason=None):
        max_s = max(forward.cigar_max_s(), reverse.cigar_max_s())
        identity_ratio = min(forward.identity_ratio, reverse.identity_ratio)
        single_best_match = (forward.alignment_scores[0] > forward.alignment_scores[1]) and \
                            (reverse.alignment_scores[0] > reverse.alignment_scores[1])
        self.conn.execute(self.table.insert(), forward_name=forward.qname, reverse_name=reverse.qname,
                          forward_sequence=forward.seq, reverse_sequence=reverse.seq, identity_ratio=identity_ratio,
                          max_s=max_s, single_best_match=single_best_match, reject_reason=reject_reason)

    def find_by_over_hang_and_identity(self, forward_file_name, reverse_file_name, max_s=None, identity=None):
        forward_file = open(forward_file_name, 'w')
        reverse_file = open(reverse_file_name, 'w')

        query = self.table.select()
        if max_s is not None:
            query = query.where(self.table.c.max_s <= max_s)
        if identity is not None:
            query = query.where(self.table.c.identity_ratio >= identity)

        results = self.conn.execute(query)
        for result in results:
            forward_file.write('{}\n{}\n'.format(result['forward_name'], result['forward_sequence']))
            reverse_file.write('{}\n{}\n'.format(result['reverse_name'], result['reverse_sequence']))

        forward_file.close()
        reverse_file.close()

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


parser = argparse.ArgumentParser(description='Generates fasta files based on qualities')
parser.add_argument('db', type=str, help='Database file location')
parser.add_argument('forward_file', type=str, help='Where forward files will be written')
parser.add_argument('reverse_file', type=str, help='Where reverse files will be written')
parser.add_argument('s_to_allow', type=float, help='Allowable query overhang')
parser.add_argument('keep_percent', type=restricted_float, help='percent to keep')

if __name__ == '__main__':
    args = parser.parse_args()
    storage  = BowtiePairStorage(args.db)
    storage.find_by_over_hang_and_identity(args.forward_file, args.reverse_file, args.s_to_allow, args.keep_percent)




#this is if you are going to be dealing with the sequences singly, not as pairs
class BowtieStorageSingle:
    def __init__(self, db_path):
        self.engine = sql.create_engine('sqlite:///' + db_path)
        self.conn = self.engine.connect()

        metadata = sql.MetaData()
        self.table = sql.Table('bowtie', metadata,
                               Column('name', Text),
                               Column('pairs_name', Text),
                               Column('sequence', Text),
                               Column('direction', Text),
                               Column('identity_ratio', Float),
                               Column('max_s', Integer),
                               Column('single_best_match', Integer),
                               Column('reject_reason', Integer)   #,
                               )

        metadata.create_all(self.engine, checkfirst=True)

    def store_pair(self, forward_entry, reverse_entry, reject_reason):
        self.store_one(forward_entry, reject_reason, 'forward', reverse_entry.qname)
        self.store_one(reverse_entry, reject_reason, 'reverse', forward_entry.qname)

    def store_one(self, entry, reject_reason=None, direction=None, pairs_name=None):
        single_best_match = entry.alignment_scores[0] > entry.alignment_scores[1]
        self.conn.execute(self.table.insert(), name=entry.qname, sequence=entry.seq,
                          pairs_name=pairs_name, direction=direction, reject_reason=reject_reason,
                          identity_ratio=entry.identity_ratio, single_best_match=single_best_match)




