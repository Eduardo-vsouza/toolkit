import os

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
from Bio import SeqIO
from venn import venn


class Intersection(object):
    def __init__(self, fasta_folder):
        self.fastaFolder = fasta_folder

        self.torfs = []
        self.gorfs = []

        self.uniqueTORFs = []
        self.uniqueGORFs = []
        self.bORFs = []

    def grab_sequences(self, tiers):
        files = os.listdir(self.fastaFolder)
        for file in files:
            if 'fasta' in file:
                print(file)
                tier = file.replace(".fasta", "").split("_")[1]
                if tier in tiers:
                    records = SeqIO.parse(f'{self.fastaFolder}/{file}', 'fasta')
                    for record in records:
                        seq = str(record.seq)
                        if 'gorfs' in file:
                            self.gorfs.append(seq)
                        elif 'torfs' in file:
                            self.torfs.append(seq)

    def plot(self):
        out = venn2([set(self.gorfs), set(self.torfs)], set_labels=('gORFs', 'tORFs'), set_colors=('#0047cb', 'green'))
        venn2_circles([set(self.gorfs), set(self.torfs)], linewidth=1, color='black', linestyle='--')
        for text in out.set_labels:
            text.set_fontsize(16)
        for text in out.subset_labels:
            text.set_fontsize(16)
        plt.show()

    def decoy_detour(self, transcriptome_database, genome_database):
        def get_seqs(db):
            seqs = []
            records = SeqIO.parse(db, 'fasta')
            for record in records:
                seqs.append(str(record.seq))
            return seqs
        t_db = get_seqs(transcriptome_database)
        g_db = get_seqs(genome_database)
        for orf in self.torfs:
            if orf in g_db:
                self.bORFs.append(orf)
            else:
                self.uniqueTORFs.append(orf)
        for orf in self.gorfs:
            if orf in t_db:
                self.bORFs.append(orf)
            else:
                self.uniqueGORFs.append(orf)
        print(self.uniqueTORFs)

    def plot_decoy_venn(self):
        subsets = (len(set(self.uniqueGORFs)), len(set(self.uniqueTORFs)), len(set(self.bORFs)))
        out = venn2(subsets=subsets, set_labels=('gORFs', 'tORFs'), set_colors=('#0047cb', 'green'))
        venn2_circles(subsets, linewidth=1, color='black', linestyle='--')
        for text in out.set_labels:
            text.set_fontsize(14)
        for text in out.subset_labels:
            text.set_fontsize(16)
        plt.show()


if __name__ == '__main__':
    # data = Intersection(fasta_folder='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/tiers_analyses/torfs_gorfs_fasta')
    # data.grab_sequences(tiers=['T1', 'T2', 'T3'])
    # # data.grab_sequences(tiers=['T3', 'T4'])
    # # data.grab_sequences(tiers=['T1', 'T2', 'T3', 'T4'])
    # # data.plot()
    # data.decoy_detour(transcriptome_database='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/transcriptome_database.fasta',
    #                   genome_database='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/genome_database.fasta')
    # data.plot_decoy_venn()

    """ re-analysis 03/04/22 after reformatting the altmeth stuff """
    # data = Intersection(fasta_folder='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/torfs_gorfs_fasta')
    # # data.grab_sequences(tiers=['T1', 'T2', 'T3'])
    # data.grab_sequences(tiers=['T1', 'T2'])
    #
    # # data.grab_sequences(tiers=['T1', 'T2', 'T3', 'T4'])
    # # data.grab_sequences(tiers=['T1', 'T2', 'T3', 'T4', 'T5'])
    # data.decoy_detour(transcriptome_database='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptome_database.fasta',
    #                   genome_database='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/genome_database.fasta')
    # # data.plot_decoy_venn()
    # # data.decoy_detour()
    # data.plot()

    """ reanalysis after applying peptide lvl fdr cutoff """
    data = Intersection(fasta_folder='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/gorfs_torfs_fasta')
    # data.grab_sequences(tiers=['T1', 'T2'])
    # data.grab_sequences(tiers=['T1', 'T2', 'T3'])
    # data.grab_sequences(tiers=['T1', 'T2', 'T3', 'T4'])
    data.grab_sequences(tiers=['T1', 'T2', 'T3', 'T4', 'T5'])
    # data.decoy_detour(transcriptome_database='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptome_database.fasta',
    #                   genome_database='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/genome_database.fasta')
    # data.plot_decoy_venn()
    # data.decoy_detour()
    data.plot()