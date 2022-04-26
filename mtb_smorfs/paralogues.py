import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class HomologyChecker:
    def __init__(self, df, tiers):
        self.df = pd.read_csv(df, sep='\t')
        self.df = self.df[self.df["Tier"].isin(tiers)]
        self.df = self.df.drop_duplicates(subset=['Final Entries'])
        self.entries = self.df["Final Entries"].tolist()
        self.sequences = self.df["chosen_sequences"].tolist()

    def check_paralogues(self):
        gorfs_paralogues = {}
        torfs_paralogues = {}
        toplot = {'Subset': [], 'Paralogues': [], 'Sequence': []}
        for i in range(len(self.entries)):
            entry = self.entries[i]
            seq = self.sequences[i]
            if 'tORF' in entry:
                if seq not in torfs_paralogues:
                    torfs_paralogues[seq] = []
                torfs_paralogues[seq].append(entry)
            elif 'gORF' in entry:
                if seq not in gorfs_paralogues:
                    gorfs_paralogues[seq] = []
                gorfs_paralogues[seq].append(entry)
        print('gORFs\n')
        for gorf in gorfs_paralogues:
            # print(gorf, len(gorfs_paralogues[gorf]))
            toplot['Subset'].append('gORF')
            toplot['Paralogues'].append(len(gorfs_paralogues[gorf]))
            toplot['Sequence'].append(gorf)
        for torf in torfs_paralogues:
            toplot['Subset'].append('tORF')
            toplot['Paralogues'].append(len(torfs_paralogues[torf]))
            toplot['Sequence'].append(torf)
            if len(torfs_paralogues[torf]) > 1:
                print(torf)
            # print(torf, len(torfs_paralogues[torf]))
        self.toPlot = toplot

    def plot(self):
        df = pd.DataFrame(data=self.toPlot)
        ax = sns.boxplot(data=df, x='Subset', y='Paralogues')
        plt.show()

if __name__ == '__main__':
    tiers = ['T1', 'T2', 'T3']
    data = HomologyChecker(df='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls', tiers=tiers)
    data.check_paralogues()
    data.plot()