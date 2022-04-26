#

import sys
import os

import pandas as pd


class SampleRenamer:
    def __init__(self, df, folder):
        self.folder = folder
        self.df = pd.read_csv(df, sep='\t')

        self.names = {}
        self.replicates = {}

    def get_names(self):
        runs = self.df["run"].tolist()
        conditions = self.df["condition"].tolist()
        for run, condition in zip(runs, conditions):
            self.names[run] = condition

    def rename_folders(self):
        files = os.listdir(self.folder)
        for file in files:
            if '.fastq' in file:
                filepath = f'{self.folder}/{file}'
                if os.path.isdir(filepath):
                    splat = file.split("_")
                    run = splat[0]
                    new_name = self.names[run]
                    print(new_name)
                    if new_name not in self.replicates:
                        self.replicates[new_name] = 0
                    self.replicates[new_name] += 1
                    final_name = f'{new_name}_{self.replicates[new_name]}'
                    cmd = f'mv {filepath} {self.folder}/{final_name}'
                    os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: .py <sample_df> <folder_with_subfolders>')
    else:
        data = SampleRenamer(df=sys.argv[1], folder=sys.argv[2])
        data.get_names()
        data.rename_folders()

