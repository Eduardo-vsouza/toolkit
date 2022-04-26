import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class IsoformPlotter(object):
    def __init__(self, tpms):
        self.tpms = pd.read_csv(tpms, sep='\t')
        self.isoforms = self.tpms["isoform_tpms"].tolist()

        self.counts = []

    def count(self):
        for i in self.isoforms:
            isos = i.split(",")
            print(isos)
            print(self.counts)
            self.counts.append(len(isos)+1)
        counts = self.counts
        counts.sort()
        counts = np.log1p(counts)
        print(counts)
        return counts



def plot(counts1, counts2):
    c = []
    group = []
    for i in counts1:
        c.append(i)
        group.append('102')
    for i in counts2:
        c.append(i)
        group.append('421')
    df = pd.DataFrame(data={'isoforms': c, 'subset': group})
    sns.histplot(data=df, x='isoforms', hue='subset', multiple='stack', kde=True)
    plt.show()


iso_421 = IsoformPlotter(tpms='/media/eduardo/DATA/Eduardo/adenovirus/nanopore/corrected_transcripts_assemblies/stringtie/421/421_tpms.txt')
counts_421 = iso_421.count()

iso_102 = IsoformPlotter(tpms='/media/eduardo/DATA/Eduardo/adenovirus/nanopore/corrected_transcripts_assemblies/stringtie/102/102_tpms.txt')
counts_102 = iso_102.count()

plot(counts_102, counts_421)


