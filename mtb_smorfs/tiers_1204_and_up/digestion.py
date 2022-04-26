import os
import sys

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class Digestion:
    def __init__(self, folder, outdir):
        self.folder = folder
        self.outdir = outdir

        self.groups = {}

    def digest(self):
        if not os.path.exists(self.outdir):
            os.system(f'mkdir {self.outdir}')
        files = os.listdir(self.folder)
        for file in files:
            cmd = f'rpg -e 42 -i {self.folder}/{file} -o {self.outdir}/{file.replace(".fasta", "")}_digested'
            os.system(cmd)

    def get_peptide_len(self):
        files = os.listdir(self.outdir)
        for file in files:
            if 'annotated' in file:
                pattern = 'Annotated'
            else:
                pattern = file.split("_")[0]
            self.groups[pattern] = []

            records = SeqIO.parse(f'{self.outdir}/{file}', 'fasta')
            for record in records:
                seq = str(record.seq)
                self.groups[pattern].append(len(seq))

    def plot(self):
        lengths = []
        groups = []
        for group in self.groups:
            for l in self.groups[group]:
                lengths.append(l)
                groups.append(group)
        df = pd.DataFrame(data={'Peptide length': lengths, 'Group': groups})
        # sns.boxplot(data=df, x='Group', y='Peptide length')
        # sns.histplot(data=df, x='Peptide length', hue='Group', multiple='dodge', common_norm = False, stat = 'density',
        #              kde=True)
        sns.displot(data=df, x='Peptide length', kind='kde', hue='Group', common_norm=False)
        # sns.lineplot(data=df, x='Peptide length', hue='Group')
        plt.xlim(6, 40)
        plt.ylim(0, 0.08)
        plt.show()

    def plot_line(self):
        lengths = []
        groups = []
        props = []

        lens = {}
        total = {}
        for group in self.groups:
            if group not in total:
                total[group] = 0
            for l in self.groups[group]:
                if group not in lens:
                    lens[group] = {}
                if l not in lens[group]:
                    lens[group][l] = 0
                lens[group][l] += 1
                total[group] += 1
        print(lens)
        for group in lens:
            for l in lens[group]:
                # print(lens[group][l])
                print('length: ', l)
                print(lens[group][l])
                props.append((lens[group][l] / total[group]) * 100)
                groups.append(group)
                lengths.append(l)

        print(lens)


        df = pd.DataFrame(data={'Peptide length': lengths, 'Group': groups, 'Proportion': props})
        sns.lineplot(data=df, x='Peptide length', y='Proportion', hue='Group')
        plt.show()




if __name__ == '__main__':
    # if sys.argv[1] == '-h' or sys.argv[1] == '--help':
    #     print("usage: .py ")
    # else:
    folder = '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/smorfs_tiers_fasta'
    # data = Digestion(folder=sys.argv[1], outdir=sys.argv[2])
    data = Digestion(folder=f'{folder}/dedupli', outdir=f'{folder}/digested')

    # data.digest()
    data.get_peptide_len()
    data.plot()
    # data.plot_line()