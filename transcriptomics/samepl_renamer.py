import os
import sys

import pandas as pd


class SampleRenamer:
    def __init__(self, df, folder):
        self.df = pd.read_csv(df, sep='\t')
        self.folder = folder

        self.groups = {}
        self.groupNumber = {}
        self.__get_groups()

    def __get_groups(self):
        samples = self.df["sample"].tolist()
        groups = self.df["group"].tolist()
        for sample, group in zip(samples, groups):
            self.groups[sample] = group
            if group not in self.groupNumber:
                self.groupNumber[group] = 0

    def rename(self):
        folders = os.listdir(self.folder)
        for folder in folders:
            sample = folder.split("_")[0]
            if sample in self.groups:
                self.groupNumber[self.groups[sample]] += 1
                cmd = f'mv {folder} {self.groups[sample]}_{self.groupNumber[self.groups[sample]]}'
                os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: <sample_df> <folder>")
    else:
        data = SampleRenamer(df=sys.argv[1], folder=sys.argv[2])
        data.rename()
