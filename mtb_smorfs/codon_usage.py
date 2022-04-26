import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class CodonUsage(object):
    def __init__(self, results):
        self.df = pd.read_csv(results, sep='\t')
        self.df = self.df.drop_duplicates(subset=['chosen_sequences'])
        self.codons = self.df["start_codon"].tolist()
        self.tiers = self.df["Tier"].tolist()
        self.entries = self.df["Final Entries"].tolist()

    def get_codons(self):
        codon_tier = {}
        i = 0
        clist = ['ATT', 'ATG', 'GTG', 'TTG']
        tiers = ['T1', 'T2', 'T3', 'T4', 'T5']
        for tier in tiers:
            codon_tier[tier] = {}
            for codon in clist:
                codon_tier[tier][codon] = 0

        for codon, tier in zip(self.codons, self.tiers):
            if codon in clist:
                if tier not in codon_tier:
                    codon_tier[tier] = {}
                if codon not in codon_tier[tier]:
                    codon_tier[tier][codon] = 0
                codon_tier[tier][codon] += 1
        print(codon_tier)
        self.codonTier = codon_tier
        return codon_tier



    def plot(self, out_folder):
        colors = sns.color_palette('pastel')[0:4]
        self.codonTier = dict(sorted(self.codonTier.items()))
        for tier in self.codonTier:
            labels = list(self.codonTier[tier].keys())
            values = list(self.codonTier[tier].values())
            plt.pie(values, labels=labels, colors=colors, autopct='%.0f%%', textprops={'fontsize': 14})
            plt.savefig(f'{out_folder}/{tier}_codon_usage.png')
            plt.clf()
            # plt.show()


class AnnoCodonUsage:
    def __init__(self, df):
        self.df = pd.read_csv(df, sep='\t')
        self.codons = self.df["start_codon"].tolist()
        print(self.df[self.df["start_codon"] == 'CTG'])

        self.codonUsage = {'ATT': 0, 'ATG': 0, 'CTG': 0, 'GTG': 0, 'ATC': 0, 'TTG': 0}

    def get_codons(self):
        for codon in self.codons:
            if codon not in self.codonUsage:
                self.codonUsage[codon] = 0
            self.codonUsage[codon] += 1
        print(self.codonUsage)

    def plot(self, out_folder):
        colors = sns.color_palette('pastel')[0:6]
        labels = list(self.codonUsage.keys())
        values = list(self.codonUsage.values())
        # plt.pie(values, labels=labels, colors=colors, autopct='%.0f%%', textprops={'fontsize': 14})
        plt.pie(values, labels=labels, colors=colors, textprops={'fontsize': 14})

        plt.savefig(f'{out_folder}/annotated_codon_usage_nopct.png')
        plt.clf()


if __name__ == '__main__':
    # data = CodonUsage(results='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/Transcriptome/Results/transcriptome_results_pre_validation_tiers_strand_coordinates_fixed_cds_starts.txt')
    # data.get_codons()
    # data.plot()

    # genome = CodonUsage(results='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/Genome/Results/genome_results_tier_with_cds.txt')
    # genome.get_codons()

    # cat = CodonUsage(results='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/analysis/genome_transcriptome_cat_results_pre_validation_codon_usage.txt')
    # cat.get_codons()
    # cat.plot(out_folder='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/analysis/codons')

    """ re-analysis 03/04/22 after reformatting the altmeth stuff """
    """ 14 04 """
    # cat = CodonUsage(results='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/cat_transcriptome_genome_new_tiers_pep_fdr.csv')
    # cat.get_codons()
    # cat.plot(out_folder='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/codons')

    anno = AnnoCodonUsage(df='/media/eduardo/gold/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info_with_rbs.txt')
    anno.get_codons()
    anno.plot('/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/codons')
