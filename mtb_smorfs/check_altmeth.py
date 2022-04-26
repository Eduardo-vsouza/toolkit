import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class AltMethInspector(object):
    def __init__(self, results, tierlist):
        self.df = pd.read_csv(results, sep='\t')
        self.df = self.df.drop_duplicates(subset=["Final Entries"])
        self.df = self.df[self.df["Tier"].isin(tierlist)]
        self.entries = self.df["Final Entries"].tolist()

        self.alternatives = {}


    def get_alternatives(self):
        for entry in self.entries:
            if 'tORF' in entry:
                splat = entry.split("_")
                name = '_'.join(splat[:3])
                original = '-'.join(name.split("-")[:3]).replace("-altmeth", "").replace("-extended", "")
                # break
                # print(name)
                # print(original, '\n')
            elif 'gORF' in entry:
                splat = entry.split("_")
                # del splat[1]
                # print(splat)
                name = '_'.join(splat[:3])

                original = '-'.join(name.split("-")).replace("-altmeth", "").replace("-extended", "")
                # print(name)
                #
                # print(entry)
                # print(original, '\n')

            if original not in self.alternatives:
                self.alternatives[original] = []
            if entry not in self.alternatives[original]:
                self.alternatives[original].append(entry)
        # print(self.alternatives)

    def plot_alternatives(self):
        subset = []
        counts = []

        for alt in self.alternatives:
            if 'tORF' in alt:
                subset.append('tORF')
            elif 'gORF' in alt:
                subset.append('gORF')
            counts.append(len(self.alternatives[alt]))
        df = pd.DataFrame(data={'Subset': subset, 'Counts': counts})
        ax = sns.histplot(data=df, x='Counts', hue='Subset', multiple='dodge', shrink=0.9)
        # ax = sns.histplot(data=df, x='Counts')

        for container in ax.containers:
            ax.bar_label(container)
        plt.xlim(0, 5)
        plt.show()


    def plot_pie(self):
        torfs_alts = 0
        gorfs_alts = 0
        torfs = 0
        gorfs = 0
        for alt in self.alternatives:
            if 'tORF' in alt:
                # subset.append('tORF')
                if len(self.alternatives[alt]) == 1:
                    torfs += 1
                else:
                    torfs_alts += len(self.alternatives[alt])
                if len(self.alternatives[alt]) >= 2 <= 3:
                    print(alt, self.alternatives[alt])
            elif 'gORF' in alt:
                if len(self.alternatives[alt]) == 1:
                    gorfs += 1
                else:
                    gorfs_alts += len(self.alternatives[alt])
        data = [gorfs, gorfs_alts, torfs, torfs_alts]
        labels = ['gORFs', 'gORFs-AltMeth', 'tORFs', 'tORFs-AltMeth']
        colors = sns.color_palette('pastel')[0:5]
        plt.pie(data, labels=labels, colors=colors, autopct='%.0f%%', textprops={'fontsize': 12})
        # plt.show()

if __name__ == '__main__':
    # data = AltMethInspector(results='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/analysis/genome_transcriptome_cat_results_pre_validation_codon_usage.txt',
    #                         tierlist=('T1', 'T2', 'T3', 'T4'))
    data = AltMethInspector(results='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/cat_genome_transcriptome_tiers_pep_fdr.xls',
                            tierlist=('T1', 'T2', 'T3', 'T4'))
    data.get_alternatives()
    # data.plot_alternatives()
    data.plot_pie()