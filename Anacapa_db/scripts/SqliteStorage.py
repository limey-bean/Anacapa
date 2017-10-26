# Written by Jesse Gomer (jessegomer@gmail.com)
# for the University of California Conservation Consortium's CALeDNA Program


import sqlalchemy as sql
from sqlalchemy import Column, Integer, Float, Text, Index
import sqlite3


class BowtieStorage:
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
                               # Index('idx_name', 'name'),
                               # Index('idx_pairs_name', 'pairs_name'),
                               # Index('idx_direction', 'direction'),
                               # Index('idx_identity_ratio', 'identity_ratio'),
                               # Index('idx_max_s', 'max_s'),
                               # Index('idx_single_best_match', 'single_best_match')
                               )

        metadata.create_all(self.engine, checkfirst=True)
        # self.conn.execute(self.table.insert(), name='abc', sequence='ATTAGAGA',
        #                   identity_ratio=0.923, max_s=100,best_match_score=120,
        #                   next_match_score=100)

    def store_pair(self, forward_entry, reverse_entry, reject_reason):
        self.store_one(forward_entry, reject_reason, 'forward', reverse_entry.qname)
        self.store_one(reverse_entry, reject_reason, 'reverse', forward_entry.qname)

    def store_one(self, entry, reject_reason=None, direction=None, pairs_name=None):
        single_best_match = entry.alignment_scores[0] <= entry.alignment_scores[1]
        self.conn.execute(self.table.insert(), name=entry.qname, sequence=entry.seq,
                          pairs_name=pairs_name, direction=direction, reject_reason=reject_reason,
                          identity_ratio=entry.identity_ratio, single_best_match=single_best_match)




