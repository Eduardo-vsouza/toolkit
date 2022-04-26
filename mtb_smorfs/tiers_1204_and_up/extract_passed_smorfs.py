import os
import sys

import pandas as pd


class ORFExtractor:
    def __init__(self, post_validation, pre_validation, pep_results_with_tiers):
        self.preValidation = pd.read_csv(pre_validation, sep='\t')
        self.postValidation = pd.read_csv(post_validation, sep='\t')
        self.results = pd.read_csv(pep_results_with_tiers, sep='\t')

        self.passedScans = []
        self.__get_peptide_fdr_scans()

        self.preValidation = self.__add_scan_info(self.preValidation)
        self.postValidation = self.__add_scan_info(self.postValidation)

    def __get_peptide_fdr_scans(self):
        """ gets the scans and mzml files that correspond to smORFs that passed the peptide FDR cutoff of 0.01. """
        files = self.results["SpecFile"].tolist()
        scans = self.results["SpecID"].tolist()
        for file, scan in zip(files, scans):
            passed = f'{file}_{scan}'
            if passed not in self.passedScans:
                self.passedScans.append(passed)

    @staticmethod
    def __add_scan_info(df):
        new_col = []
        files = df["SpecFile"].tolist()
        scans = df["SpecID"].tolist()
        for file, scan in zip(files, scans):
            passed = f'{file}_{scan}'
            new_col.append(passed)
        df.insert(0, "full_scan", new_col)
        return df

    def refilter_results(self, outdir):
        if not os.path.exists(outdir):
            os.system(f'mkdir {outdir}')
        self.preValidation = self.preValidation[self.preValidation["full_scan"].isin(self.passedScans)]
        self.postValidation = self.postValidation[self.postValidation["full_scan"].isin(self.passedScans)]

        self.preValidation.to_csv(f'{outdir}/pre_validation_with_pep_fdr.txt', sep='\t', index=False)
        self.postValidation.to_csv(f'{outdir}/post_validation_with_pep_fdr.txt', sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: .py <pep_results_with_tiers> <post_validation> <pre_validation> <outdir>')
    else:
        data = ORFExtractor(pep_results_with_tiers=sys.argv[1], post_validation=sys.argv[2], pre_validation=sys.argv[3])
        data.refilter_results(outdir=sys.argv[4])