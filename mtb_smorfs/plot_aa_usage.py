import os

import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class AAComposition(object):
    def __init__(self, fasta_folder, reference_proteome):
        """ plots stacked histograms showing the AA composition of annotated smorfs, novel tORFs, novel gORFs
        and separates them by tier"""
        self.folder = fasta_folder
        self.proteome = reference_proteome

        self.tiers = {}
        self.annotatedAAs = {}

        self.lenByTiers = {'Annotated': []}

    def get_aa_by_tier(self):
        tiers = {}
        files = os.listdir(self.folder)
        checker = {}
        for file in files:
            if '_T' in file:
                tier = file.replace(".fasta", "").split("_")[1]
                if tier not in self.tiers:
                    self.tiers[tier] = {}
                filename = f'{self.folder}/{file}'
                records = SeqIO.parse(filename, 'fasta')
                for record in records:
                    seq = str(record.seq)
                    if tier not in checker:
                        checker[tier] = []

                    if seq not in checker[tier]:
                        checker[tier].append(seq)
                        self.__get_aa(seq, tier)
                        self.__get_lens(seq, tier)
        # self.tiers = tiers

    def __get_lens(self, seq, tier):
        if tier not in self.lenByTiers:
            self.lenByTiers[tier] = []
        self.lenByTiers[tier].append(len(seq))


    def __get_aa(self, seq, tier):
        for aa in seq:
            if aa not in self.tiers[tier]:
                self.tiers[tier][aa] = 0
            self.tiers[tier][aa] += 1

    def to_log(self):
        for tier in self.tiers:
            for aa in self.tiers[tier]:
                self.tiers[tier][aa] = np.log(self.tiers[tier][aa])

    def to_proportion(self):
        for tier in self.tiers:
            total = 0
            for aa in self.tiers[tier]:
                total += self.tiers[tier][aa]
            for aa in self.tiers[tier]:
                self.tiers[tier][aa] = (self.tiers[tier][aa]/total)*100

        total = 0
        for aa in self.annotatedAAs:
            total += self.annotatedAAs[aa]
        for aa in self.annotatedAAs:
            self.annotatedAAs[aa] = (self.annotatedAAs[aa]/total)*100


    def get_annotated(self):
        records = SeqIO.parse(self.proteome, 'fasta')
        for record in records:
            seq = str(record.seq)
            for aa in seq:
                if aa != '+':
                    if aa not in self.annotatedAAs:
                        self.annotatedAAs[aa] = 0
                    self.annotatedAAs[aa] += 1
            self.lenByTiers['Annotated'].append(len(seq))

    def plot(self):
        aas = []
        tiers = []
        counts = []
        for tier in self.tiers:
            for aa in self.tiers[tier]:
                aas.append(aa)
                tiers.append(tier)
                counts.append(self.tiers[tier][aa])
        for anno in self.annotatedAAs:
            aas.append(anno)
            tiers.append('Annotated')
            counts.append(self.annotatedAAs[anno])

        df = pd.DataFrame(data={'Amino acid': aas, 'Tier': tiers, 'Proportion': counts})
        sns.color_palette("viridis", as_cmap=True)

        sns.barplot(data=df, x="Amino acid", y='Proportion', hue='Tier', palette='colorblind')
        plt.show()

    def plot_length(self):
        lens = []
        tiers = []
        for tier in self.lenByTiers:
            for l in self.lenByTiers[tier]:
                lens.append(l)
                tiers.append(tier)
        y = 'Microprotein length (aa)'
        x = 'Microprotein group'
        df = pd.DataFrame(data={y: lens, x: tiers})
        print(df[y])

        sns.boxplot(data=df, x=x, y=y)
        plt.show()

if __name__ == '__main__':
    # data = AAComposition(fasta_folder='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/tiers_analyses/torfs_gorfs_fasta',
    #                      reference_proteome='/media/eduardo/New Volume/Eduardo/mtb_annotation_files/mtb_reference_smorfome.fasta')
    # data.get_aa_by_tier()
    # data.get_annotated()
    # # data.to_log()
    # data.to_proportion()
    # # data.plot()
    # data.plot_length()

    """ re-analysis 03/04/22 after reformatting the altmeth stuff 
    '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls"""
    # data = AAComposition(fasta_folder='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/torfs_gorfs_fasta',
    #                      reference_proteome='/media/eduardo/gold/Eduardo/mtb_annotation_files/mtb_reference_smorfome.fasta')
    # data.get_aa_by_tier()
    # data.get_annotated()
    # data.to_proportion()
    # # data.plot_length()
    # data.plot()

    """ reanalysis 12/04/22 after applying peptide fdr cutoff """
    data = AAComposition(fasta_folder='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/gorfs_torfs_fasta',
                         reference_proteome='/media/eduardo/gold/Eduardo/mtb_annotation_files/mtb_reference_smorfome.fasta')
    data.get_aa_by_tier()
    data.get_annotated()
    data.to_proportion()
    # data.plot_length()
    data.plot()
