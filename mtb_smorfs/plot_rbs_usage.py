import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


class RBSUsage(object):
    def __init__(self, anno_df, cat_results_04_genome_transcriptome):
        self.annotated = pd.read_csv(anno_df, sep='\t')
        self.novel = pd.read_csv(cat_results_04_genome_transcriptome, sep='\t')
        self.novel = self.novel.drop_duplicates(subset=["chosen_sequences"])
        self.novel = self.novel[self.novel["Tier"] != "Tier"]
        self.novel = self.novel[self.novel["Free Energy"] != "Free Energy"]
        print(self.novel["Free Energy"])
        self.rbs = {'Annotated': [], 'T1': [], 'T2': [], 'T3': [], 'T4': [], 'T5': []}

        self.rbsPresence = {'Annotated': {}, 'T1': {}, 'T2': {}, 'T3': {}, 'T4': {}, 'T5': {}}

    def get_rbs_novel(self):
        tier = self.novel["Tier"].tolist()
        shine = self.novel["Shine Dalgarno"].tolist()
        energies = self.novel["Free Energy"].tolist()

        for i in range(len(tier)):
            if tier[i] not in self.rbs:
                self.rbs[tier[i]] = []
            self.rbs[tier[i]].append(float(energies[i]))
            if shine[i] not in self.rbsPresence[tier[i]]:
                self.rbsPresence[tier[i]][shine[i]] = 0
            self.rbsPresence[tier[i]][shine[i]] += 1

    def get_rbs_annotated(self):
        energies = self.annotated["Free Energy"].tolist()
        shine = self.annotated["Shine Dalgarno"].tolist()
        for i in range(len(energies)):
            self.rbs['Annotated'].append(float(energies[i]))
            if shine[i] not in self.rbsPresence['Annotated']:
                self.rbsPresence['Annotated'][shine[i]] = 0
            self.rbsPresence['Annotated'][shine[i]] += 1

    def plot_energies(self):
        energies = []
        subsets = []
        for tier in self.rbs:
            for energy in self.rbs[tier]:
                subsets.append(tier)
                energies.append(energy)
        df = pd.DataFrame(data={'Free Energy': energies,
                                'Subset': subsets})
        print(df["Free Energy"])
        sns.boxplot(data=df, x='Subset', y='Free Energy')
        plt.show()

    def shine_proportions(self):
        proportions = {}
        for tier in self.rbsPresence:
            total = 0
            for shine in self.rbsPresence[tier]:
                total += self.rbsPresence[tier][shine]
            for shine in self.rbsPresence[tier]:
                proportions[tier] = self.rbsPresence[tier]
                # proportions[]
                proportions[tier][shine] = (self.rbsPresence[tier][shine] / total) * 100
                # print(self.rbsPresence[tier])
        self.proportions = proportions

    def plot_shine(self):
        shines = []
        subsets = []
        presence = []
        for tier in self.rbsPresence:
            for shine in self.rbsPresence[tier]:
                shines.append(shine)

                subsets.append(tier)
                presence.append(self.rbsPresence[tier][shine])
        # presence = np.log1p(presence)
        fixed_shines = []
        for i in shines:
            if i == 'Leaderless':
                fixed = 'SD absent'
            else:
                fixed = i
            fixed_shines.append(fixed)
        df = pd.DataFrame(data={'Shine Dalgarno': fixed_shines, 'Subset': subsets, 'RBS counts': presence})
        sns.barplot(data=df, x='Subset', y='RBS counts', hue='Shine Dalgarno', edgecolor='black')
        plt.show()

    def plot_shine_proportions(self):
        self.shine_proportions()
        shines = []
        subsets = []
        presence = []
        for tier in self.proportions:
            for shine in self.proportions[tier]:
                shines.append(shine)

                subsets.append(tier)

                presence.append(self.proportions[tier][shine])
        # presence = np.log1p(presence)
        fixed_shines = []
        for i in shines:
            if i == 'Leaderless':
                fixed = 'Absent.'
            else:
                if i == 'SD absent':
                    fixed = 'Absent.'
                else:
                    fixed = i
            fixed_shines.append(fixed)
        df = pd.DataFrame(data={'Shine Dalgarno': fixed_shines, 'Subset': subsets, 'RBS (%)': presence})
        sns.barplot(data=df, x='Subset', y='RBS (%)', hue='Shine Dalgarno', edgecolor='black')
        # plt.tight_layout()
        plt.ylim(0, 80)
        plt.show()


if __name__ == '__main__':
    # data = RBSUsage(cat_results_04_genome_transcriptome='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/analysis/genome_transcriptome_cat_results_pre_validation_codon_usage.txt',
    #                 anno_df='/media/eduardo/New Volume/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info_with_rbs.txt')
    # data.get_rbs_novel()
    # data.get_rbs_annotated()
    # # data.plot()
    # # data.plot_shine()
    # data.plot_shine_proportions()

    """ re-analysis 03/04/22 after reformatting the altmeth stuff 
    '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls"""
    # data = RBSUsage(cat_results_04_genome_transcriptome='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls',
    #                 anno_df='/media/eduardo/gold/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info_with_rbs.txt')
    # data.get_rbs_novel()
    # data.get_rbs_annotated()
    # data.plot_shine_proportions()
    # data.plot_shine()
    # data.plot_energies()

    """ reanalysis 14/03/22 """
    data = RBSUsage(cat_results_04_genome_transcriptome='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/cat_transcriptome_genome_new_tiers_pep_fdr.csv',
                    anno_df='/media/eduardo/gold/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info_with_rbs.txt')
    data.get_rbs_novel()
    data.get_rbs_annotated()
    # data.plot_shine_proportions()
    data.plot_shine()
    data.plot_energies()
