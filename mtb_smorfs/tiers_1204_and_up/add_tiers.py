import os
import sys

import pandas as pd


class TierInclusion:
    def __init__(self, results, tier):
        self.results = pd.read_csv(results, sep='\t')
        self.tierDF = pd.read_csv(tier, sep='\t')

        self.tiers = {}
        self.__get_seq_tier()

    def __get_seq_tier(self):
        tiers = self.tierDF["tiers"].tolist()
        seqs = self.tierDF["orfs"].tolist()
        for tier, seq in zip(tiers, seqs):
            self.tiers[seq] = tier

    def refilter(self, output):
        tiers = []  # tier in the results data frame. to be inserted
        seqs = self.results["chosen_sequences"].tolist()
        for seq in seqs:
            tiers.append(self.tiers[seq])
        self.results = self.results.drop(columns=['Tier'])
        self.results.insert(4, "Tier", tiers)
        self.results.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: <results> <tier_df> <output>')
    else:
        data = TierInclusion(results=sys.argv[1], tier=sys.argv[2])
        data.refilter(output=sys.argv[3])