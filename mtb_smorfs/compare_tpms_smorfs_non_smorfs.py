import os
import sys

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


class TPMAnalysis(object):
    def __init__(self, folder, transcriptome_results):
        """

        :param folder: folder containing subfolders where the quant.sf files from Salmon output are stored
        :param transcriptome_results: yup
        """
        self.folder = folder
        self.df = pd.read_csv(transcriptome_results, sep='\t')

        self.scts = []
        self.normalTranscripts = []

        self.sctsTPMs = {}
        self.normalTranscriptsTPMs = {}

        self.meansSCTs = {}
        self.meansNormal = {}

    def identify_scts(self, tiers):
        """ identifies smORF-containing transcripts (SCTs) """
        self.df = self.df[self.df["Tier"].isin(tiers)]
        entries = self.df["Final Entries"].tolist()
        for entry in entries:
            splat = entry.split("_")
            gene = splat[1].split("-")
            transcript = '-'.join(gene[:2])
            if transcript not in self.scts:
                self.scts.append(transcript)

    def get_tpms(self, pattern):
        subdirs = os.listdir(self.folder)
        for folder in subdirs:
            if os.path.isdir(f'{self.folder}/{folder}'):
                if pattern in folder:
                    quant = f'{self.folder}/{folder}/quant.sf'
                    df = pd.read_csv(quant, sep='\t')
                    transcripts = df["Name"].tolist()
                    tpms = df["TPM"].tolist()
                    replicate = folder.split(".")[0]
                    for transcript, tpm in zip(transcripts, tpms):
                        if transcript in self.scts:
                            if transcript not in self.sctsTPMs:
                                self.sctsTPMs[transcript] = {}
                            self.sctsTPMs[transcript][replicate] = int(tpm)
                        else:
                            if transcript not in self.normalTranscriptsTPMs:
                                self.normalTranscriptsTPMs[transcript] = {}
                            self.normalTranscriptsTPMs[transcript][replicate] = int(tpm)
        # print(self.sctsTPMs)
        # print(self.normalTranscriptsTPMs)

    def get_means(self):
        for transcript in self.sctsTPMs:
            reps = []
            for rep in self.sctsTPMs[transcript]:
                reps.append(self.sctsTPMs[transcript][rep])
            self.meansSCTs[transcript] = np.mean(reps)
        for transcript in self.normalTranscriptsTPMs:
            reps = []
            for rep in self.normalTranscriptsTPMs[transcript]:
                reps.append(self.normalTranscriptsTPMs[transcript][rep])
            self.meansNormal[transcript] = np.mean(reps)

    def plot(self, output):
        data = {'Transcript subset': [], 'TPMs': []}

        for transcript in self.meansSCTs:
            data['Transcript subset'].append('SCT')
            data['TPMs'].append(self.meansSCTs[transcript])
        for transcript in self.meansNormal:
            data['Transcript subset'].append('Annotated')
            data['TPMs'].append(self.meansNormal[transcript])
        df = pd.DataFrame(data=data)
        df.to_csv(output, sep='\t', index=False)
        sns.barplot(data=df, x='Transcript subset', y='TPMs', edgecolor='black')
        # sns.boxplot(data=df, x='Transcript subset', y='TPMs')

        plt.ylabel('Mean TPM values')
        # plt.tight_layout()
        plt.show()



if __name__ == '__main__':
    # data = TPMAnalysis(folder='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/reads/e-mtab-1616/quantification',
    #                    transcriptome_results='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/Transcriptome/Results/transcriptome_results_pre_validation_tiers_strand_coordinates_fixed_cds_starts.txt')

    # 13/04/22
    data = TPMAnalysis(
        folder='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/quantification/mtab-1616',
        transcriptome_results='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/Transcriptome/Results/reanalysis_pep_fdr/transcriptome_results_with_tiers_pep_fdr_tiers_redone.xls')

    data.identify_scts(tiers=('T1', 'T2', 'T3'))
    data.get_tpms(pattern='control')
    data.get_means()
    data.plot(output='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/quantification/control_TPM_scts_and_normal_rnas.xls')
