#!/usr/bin/env python3

import os
import sys

import pandas as pd


class HCFilter(object):
    def __init__(self, reformatted, hc_results):
        self.reformatted = pd.read_csv(reformatted, sep='\t')
        self.hcResults = pd.read_csv(hc_results, sep='\t')

    def filter_scans(self, output):
        files = self.hcResults["SpecFile"].tolist()
        specs = self.hcResults["SpecID"].tolist()
        self.reformatted = self.reformatted[self.reformatted["SpecFile"].isin(files)]
        self.reformatted = self.reformatted[self.reformatted["SpecID"].isin(specs)]
        self.reformatted.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: .py <reformatted_results> <hc_results> <output>')
    else:
        data = HCFilter(reformatted=sys.argv[1], hc_results=sys.argv[2])
        data.filter_scans(output=sys.argv[3])
