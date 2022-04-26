import os
import sys
import pandas as pd


class SalmonCountRenamer(object):
    def __init__(self, count_file):
        self.df = pd.read_csv(count_file, sep='\t')

    def check_alts(self):
        names = self.df["Name"].tolist()
        for name in names:
            if