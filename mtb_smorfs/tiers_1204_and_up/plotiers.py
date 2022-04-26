import os
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class TierPlotting:
    def __init__(self, transcriptome, genome):
        self.genome = pd.read_csv(genome, sep='\t')
        self.transcriptome = pd.read_csv(transcriptome, sep='\t')

        self.genome = self.genome.drop_duplicates(subset=["chosen_sequences"])
        self.transcriptome = self.transcriptome.drop_duplicates(subset=["chosen_sequences"])

        self.df = pd.DataFrame()

    def count_tiers(self):
        tiers = []
        subsets = []
        counts = []

        def add(df, subset):
            tier = df["Tier"].tolist()
            tier_counts = {t: 0 for t in set(tier)}

            for t in tier:
                tier_counts[t] += 1
            for t in tier_counts:
                tiers.append(t)
                subsets.append(subset)
                counts.append(tier_counts[t])

        add(self.transcriptome, "Transcriptome")
        add(self.genome, "Genome")
        data = {'Tier': tiers, 'Subset': subsets, 'smORFs': counts}

        df = pd.DataFrame(data=data)

        colors = ['#5b96d6', '#5dbf70']
        ax = sns.barplot(x="Tier", y="smORFs", hue="Subset", data=df, edgecolor='black',
                         order=['T1', 'T2', 'T3', 'T4', 'T5'],
                         palette=colors[::-1])
        for container in ax.containers:
            ax.bar_label(container)
        plt.show()

if __name__ == '__main__':
    data = TierPlotting(genome=sys.argv[1], transcriptome=sys.argv[2])
    data.count_tiers()
