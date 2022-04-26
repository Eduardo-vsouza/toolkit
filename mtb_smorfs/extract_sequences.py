#!/usr/bin/env python3

import sys

import pandas as pd


class SequenceExtractor(object):
    def __init__(self, df):
        self.df = pd.read_csv(df, sep='\t')
        self.df = self.df[self.df["chosen_sequences"].str.len() <= 100]

    def to_fasta(self, output):
        sequences = self.df["chosen_sequences"].tolist()
        entries = self.df["Final Entries"].tolist()
        tiers = self.df["Tier"].tolist()
        fastas = {}
        checker = []
        for tier in tiers:
            if tier not in fastas:
                fastas[tier] = []
        for i in range(len(entries)):
            entry, seq = entries[i], sequences[i]
            tier = tiers[i]
            if entry not in checker:
                checker.append(entry)
                fastas[tier].append(f'>{entry}\n{seq}\n')
        for tier in fastas:
            with open(f'{output}_{tier}.fasta', 'w') as out:
                out.writelines(fastas[tier])

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: extract_sequences.py <results_df> <output_fasta_file>\n'
              'Do not include the ".fasta" at the end of the output file name')
    else:
        data = SequenceExtractor(df=sys.argv[1])
        data.to_fasta(output=sys.argv[2])