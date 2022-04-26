import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, anderson, normaltest, mannwhitneyu, wilcoxon, kruskal
import numpy as np


class TPMPlotter:
    def __init__(self, control, treated):
        self.control = pd.read_csv(control, sep='\t')
        self.treated = pd.read_csv(treated, sep='\t')

        self.mergedDataFrame = self.__merge()

    def __merge(self):
        rna_subsets = []
        tpms = []
        groups = []

        def add(df, group):
            sub = df["Transcript subset"].tolist()
            tpm_df = df["TPMs"].tolist()
            for i in range(len(sub)):
                rna_subsets.append(sub[i])
                tpms.append(np.log1p(tpm_df[i]))
                # tpms.append(tpm_df[i])

                groups.append(group)

        add(self.control, 'control')
        add(self.treated, 'starvation')
        df = pd.DataFrame(data={'Transcript subset': rna_subsets, 'TPMs': tpms, 'Group': groups})
        return df

    def check_normality(self):
        stat, p = shapiro(self.mergedDataFrame["TPMs"].tolist())
        print('Statistics=%.3f, p=%.3f' % (stat, p))
        alpha = 0.05
        if p > alpha:
            print('Sample looks Gaussian (fail to reject H0)')
        else:
            print('Sample does not look Gaussian (reject H0)')

    def check_normality_dagostino(self):
        stat, p = normaltest(self.mergedDataFrame["TPMs"].tolist())
        print('Statistics=%.3f, p=%.3f' % (stat, p))
        # interpret
        alpha = 0.05
        if p > alpha:
            print('Sample looks Gaussian (fail to reject H0)')
        else:
            print('Sample does not look Gaussian (reject H0)')

    def check_normality_anderson(self):
        result = anderson(self.mergedDataFrame["TPMs"].tolist())
        plt.hist(self.mergedDataFrame["TPMs"].tolist())
        plt.show()
        print('Statistic: %.3f' % result.statistic)
        p = 0
        for i in range(len(result.critical_values)):
            sl, cv = result.significance_level[i], result.critical_values[i]
            if result.statistic < result.critical_values[i]:
                print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
            else:
                print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))

    def mann_whitney(self):
        print(self.mergedDataFrame.columns)
        self.mergedDataFrame = self.mergedDataFrame[self.mergedDataFrame["Group"] == 'starvation']
        scts = self.mergedDataFrame[self.mergedDataFrame["Transcript subset"] == 'SCT']
        # scts = np.log1p(scts["TPMs"].tolist())
        scts = scts["TPMs"].tolist()

        normal = self.mergedDataFrame[self.mergedDataFrame["Transcript subset"] == 'Annotated']
        # normal = np.log1p(normal["TPMs"].tolist())
        normal = normal["TPMs"].tolist()

        # stat, p = mannwhitneyu(scts, normal)
        stat, p = kruskal(scts, normal)

        print('Statistics=%.10f, p=%.10f' % (stat, p))
        # interpret
        alpha = 0.05
        if p > alpha:
            print('Same distribution (fail to reject H0)')
        else:
            print('Different distribution (reject H0)')

    def plot(self):
        sns.boxplot(data=self.mergedDataFrame, x='Group', y='TPMs', hue='Transcript subset')
        plt.show()


if __name__ == '__main__':
    folder = '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/quantification/mtab-1616'
    data = TPMPlotter(control=f'{folder}/control_TPM_scts_and_normal_rnas.xls',
                      treated=f'{folder}/starvation_TPM_scts_and_normal_rnas.xls')
    # data.check_normality()
    # data.check_normality_anderson()
    # data.check_normality_dagostino()
    data.mann_whitney()
    # data.plot()