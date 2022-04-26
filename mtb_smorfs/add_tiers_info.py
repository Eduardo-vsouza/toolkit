#!/usr/bin/env python3

import sys

import pandas as pd


class TierInformation(object):
    def __init__(self, pre_validation_results, tiers):
        self.preValidation = pd.read_csv(pre_validation_results, sep='\t')
        self.tiers = pd.read_csv(tiers, sep='\t')

        self.tierInfo = {}
        self.__get_tiers()

    def __get_tiers(self):
        orfs = self.tiers["orfs"].tolist()
        tiers = self.tiers["tiers"].tolist()

        for orf, tier in zip(orfs, tiers):
            self.tierInfo[orf] = tier

    def __add_info(self, df):
        tier_col = []
        entries = df["Final Entries"].tolist()
        for entry in entries:
            if entry in self.tierInfo:
                tier = self.tierInfo[entry]
                tier_col.append(tier)
            else:
                tier_col.append('no tier')
        df.insert(4, "Tier", tier_col)
        return df

    def add_tiers(self, output):
        df = self.__add_info(self.preValidation)
        df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: add_tiers_info <pre_validation_results> <tiers_data_frame> <output>')
    else:
        data = TierInformation(pre_validation_results=sys.argv[1], tiers=sys.argv[2])
        data.add_tiers(output=sys.argv[3])
