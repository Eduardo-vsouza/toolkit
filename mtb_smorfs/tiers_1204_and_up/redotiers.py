import sys

import pandas as pd


class TierClassification:
    def __init__(self, results, post_validation):
        self.postValidation = pd.read_csv(post_validation, sep='\t')
        self.results = pd.read_csv(results, sep='\t')

        self.scansHC = []
        self.HCsmorfs = {}

        self.__get_hc_scans()
        self.__define_hc()

        self.tiers = {}

    def __get_hc_scans(self):
        files = self.postValidation["SpecFile"].tolist()
        scans = self.postValidation["SpecID"].tolist()
        for file, scan in zip(files, scans):
            self.scansHC.append(f'{file}_{scan}')

    def __define_hc(self):
        seqs = self.results["chosen_sequences"].tolist()
        files = self.results["SpecFile"].tolist()
        scans = self.results["SpecID"].tolist()
        preds = {}
        for seq in seqs:
            df = self.results[self.results["chosen_sequences"] == seq]
            df_seqs = df["chosen_sequences"].tolist()
            df_files = df["SpecFile"].tolist()
            df_scans = df["SpecID"].tolist()
            for i in range(len(df_scans)):
                full = f'{df_files[i]}_{df_scans[i]}'
                if full in self.scansHC:
                    preds[df_seqs[i]] = 'hc'
        for seq in seqs:
            if seq not in preds:
                preds[seq] = 'lc'
        self.HCsmorfs = preds

    def classify_results(self):
        predictions = []
        full_scans = []
        files = self.results["SpecFile"].tolist()
        scans = self.results["SpecID"].tolist()
        seqs = self.results["chosen_sequences"].tolist()
        for i in range(len(scans)):
            full = f'{files[i]}_{scans[i]}'
            predictions.append(self.HCsmorfs[seqs[i]])
            # if full in self.scansHC:
            #     predictions.append('hc')
            # else:
            #     predictions.append('lc')
            full_scans.append(full)
        self.results.insert(5, "prediction", predictions)
        self.results.insert(2, "full_scan", full_scans)

    def define_tiers(self):
        seqs = self.results["chosen_sequences"].tolist()
        scans = self.results["SpecFile"].tolist()
        predictions = self.results["prediction"].tolist()

        replicates = {}
        preds = {}

        for i in range(len(seqs)):
            if seqs[i] not in replicates:
                replicates[seqs[i]] = []
                preds[seqs[i]] = predictions[i]
            if scans[i] not in replicates[seqs[i]]:
                replicates[seqs[i]].append(scans[i])

        for seq in replicates:
            if preds[seq] == 'hc':
                if len(replicates[seq]) > 3:
                    self.tiers[seq] = 'T1'
                else:
                    self.tiers[seq] = 'T2'
            else:
                if len(replicates[seq]) > 3:
                    self.tiers[seq] = 'T3'
                elif len(replicates[seq]) == 3:
                    self.tiers[seq] = 'T4'
                else:
                    self.tiers[seq] = 'T5'

    def add_tiers(self, output):
        seqs = self.results["chosen_sequences"].tolist()
        tiers = []
        for seq in seqs:
            tier = self.tiers[seq]
            tiers.append(tier)
        self.results = self.results.drop(columns='Tier')
        self.results.insert(4, "Tier", tiers)
        self.results.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        ...
    else:
        data = TierClassification(post_validation=sys.argv[1], results=sys.argv[2])
        data.classify_results()
        data.define_tiers()
        data.add_tiers(output=sys.argv[3])




