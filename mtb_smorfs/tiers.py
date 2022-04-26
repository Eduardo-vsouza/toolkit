#!/usr/bin/env python3
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class ORFTiers(object):
    def __init__(self, pre_validation, post_validation):
        self.all = pd.read_csv(pre_validation, sep='\t')
        self.all = self.all[self.all["chosen_sequences"].str.len() <= 100]
        self.hc = pd.read_csv(post_validation, sep='\t')
        self.hc = self.hc[self.hc["chosen_sequences"].str.len() <= 100]

        self.allReplicates = self.__count_replicates(self.all)

        # these two are the ones used to define tiers
        self.hcReplicates = self.__count_replicates(self.hc)
        self.lcReplicates = self.__check_overlap()

    @staticmethod
    def __count_replicates(results):
        replicates = {}
        entries = results["Final Entries"].tolist()
        files = results["SpecFile"].tolist()
        ids = results["SpecID"].tolist()
        for i in range(len(entries)):
            entry = entries[i]
            if entry not in replicates:
                replicates[entry] = []
            scan = f'{files[i]}: {ids[i]}'
            if scan not in replicates[entry]:
                replicates[entry].append(scan)
        num_reps = {}
        for rep in replicates:
            num_reps[rep] = len(replicates[rep])
        return num_reps

    def __check_overlap(self):
        """ removes smORFS from the 'all' subset if they appear in the 'HC' subset. This is important to avoid
        classifying a HC as LC. """
        checked = {}
        for orf in self.allReplicates:
            if orf not in self.hcReplicates:
                checked[orf] = self.allReplicates[orf]
        return checked

    def get_tiers(self):

        tier = {'T1': [], 'T2': [], 'T3': [], 'T4': [], 'T5': []}
        for orf in self.hcReplicates:
            if self.hcReplicates[orf] >= 2:
                tier['T1'].append(orf)
            else:
                tier['T2'].append(orf)
        for orf in self.lcReplicates:
            if self.lcReplicates[orf] > 3:
                tier['T3'].append(orf)
            elif self.lcReplicates[orf] == 3:
                tier['T4'].append(orf)
            else:
                tier['T5'].append(orf)
        self.tiers = tier
        return tier

    def save_tiers(self, output):
        orfs = []
        tiers = []
        for tier in self.tiers:
            for orf in self.tiers[tier]:
                orfs.append(orf)
                tiers.append(tier)
        df = pd.DataFrame(data={'orfs': orfs, 'tiers': tiers})
        df.to_csv(output, sep='\t', index=False)



if __name__ == '__main__':

    if sys.argv[1] == '-h':
        print('Adds tier information to the data and creates a concatenated results, including all possible tiers '
              'for those two files together. \n'
              'usage: tiers.py <post_validation_results> <pre_validation results> <output_with_tiers>')
    else:
        data = ORFTiers(post_validation=sys.argv[1], pre_validation=sys.argv[2])

        tiers = data.get_tiers()
        data.save_tiers(sys.argv[3])
        for tier in tiers:
            print(tier, len(tiers[tier]))
        print('\n')


    # folder = '/media/eduardo/DATA/Eduardo/smorfs_mtb_090222/Transcriptome/Results'
    # data = ORFTiers(post_validation=f'{folder}/transcriptome_post_validation_results.txt',
    #                 pre_validation=f'{folder}/../post_perc/transcriptome_results_04_altmeth_coordinates.txt')
    #
    # tiers = data.get_tiers()
    # data.save_tiers(f'{folder}/tiers.xls')
    # for tier in tiers:
    #     print(tier, len(tiers[tier]))
    # print('\n')


    # folder = '/media/eduardo/DATA/Eduardo/smorfs_mtb_090222/Genome/Results'
    # data = ORFTiers(post_validation=f'{folder}/genome_post_validation_results.txt',
    #                 pre_validation=f'{folder}/../post_perc/genome_results_04.txt')
    # g_tiers = data.get_tiers()
    # data.save_tiers(f'{folder}/tiers.xls')
    # for tier in g_tiers:
    #     print(tier, len(g_tiers[tier]))