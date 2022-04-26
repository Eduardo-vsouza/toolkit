import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class GC(object):
    def __init__(self, anno_df, cat_results_04_genome_transcriptome):
        self.annotated = pd.read_csv(anno_df, sep='\t')
        self.novel = pd.read_csv(cat_results_04_genome_transcriptome, sep='\t')
        self.novel = self.novel[self.novel["Tier"] != "Tier"]
        self.novel = self.novel.drop_duplicates(subset=["chosen_sequences"])
        self.novel = self.novel[self.novel["Free Energy"] != "Free Energy"]
        self.gc = {'Annotated': [], 'T1': [], 'T2': [], 'T3': [], 'T4': [], 'T5': []}

    def get_percentages(self):
        orfs = self.novel["orf_sequence_nucleotides"].tolist()
        tiers = self.novel["Tier"].tolist()
        for orf, tier in zip(orfs, tiers):
            if type(orf) == str:
                gc = self.__count_gc(orf)
                self.gc[tier].append(gc)

        annorfs = self.annotated["cds"].tolist()
        for orf in annorfs:
            gc = self.__count_gc(orf)
            self.gc["Annotated"].append(gc)

    @staticmethod
    def __count_gc(seq):
        gc = 0
        for n in seq:
            if n == 'G' or n == 'C':
                gc += 1
        return (gc/len(seq)) * 100

    def plot(self):
        tiers = []
        gcs = []
        for tier in self.gc:
            for gc in self.gc[tier]:
                tiers.append(tier)
                gcs.append(gc)
        df = pd.DataFrame(data={'Subset': tiers, 'GC (%)': gcs})
        sns.boxplot(data=df, x='Subset', y='GC (%)')
        # plt.ylim(0, 100)
        plt.show()


if __name__ == '__main__':
    # data = GC(cat_results_04_genome_transcriptome='/media/eduardo/New Volume/Eduardo/smorfs_mtb_090222/analysis/genome_transcriptome_cat_results_pre_validation_codon_usage.txt',
    #           anno_df='/media/eduardo/New Volume/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info_with_rbs.txt')
    # data.get_percentages()
    # data.plot()

    """ re-analysis 03/04/22 after reformatting the altmeth stuff 
    '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls"""

    # data = GC(cat_results_04_genome_transcriptome='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls',
    #           anno_df='/media/eduardo/gold/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info_with_rbs.txt')
    # data.get_percentages()
    # data.plot()
    #
    """ reanalsis 12/04/22 after applying peptide fdr """
    data = GC(cat_results_04_genome_transcriptome='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/cat_transcriptome_genome_new_tiers_pep_fdr.csv',
              anno_df='/media/eduardo/gold/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info_with_rbs.txt')
    data.get_percentages()
    data.plot()