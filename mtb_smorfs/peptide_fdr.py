import sys

import pandas as pd


class PeptideFDR:
    def __init__(self, genome_pep_cat, transcriptome_pep_cat):
        self.genomePeptides = pd.read_csv(genome_pep_cat, sep='\t')
        self.genomePeptides = self.genomePeptides[self.genomePeptides["q-value"] != "q-value"]
        self.genomePeptides["q-value"] = pd.to_numeric(self.genomePeptides["q-value"], downcast='float')

        self.transcriptomePeptides = pd.read_csv(transcriptome_pep_cat, sep='\t')
        self.transcriptomePeptides = self.transcriptomePeptides[self.transcriptomePeptides["q-value"] != "q-value"]
        self.transcriptomePeptides["q-value"] = pd.to_numeric(self.transcriptomePeptides["q-value"], downcast='float')

        self.passedPeptides = []
        # self.__filter_peptides(self.genomePeptides)
        # self.__filter_peptides(self.transcriptomePeptides)

    def __filter_peptides(self, df):
        df = df[df["q-value"] <= 0.01]
        # print(df)
        # print(df.columns)
        peptides = df["peptide"].tolist()
        fixed_peptides = []
        for pep in peptides:
            pepf = pep.replace(".", "").replace("[UNIMOD:", "").replace("]", "").replace("-", "")
            fixed = ''.join([i for i in pepf if not i.isdigit()])
            fixed_peptides.append(fixed)
        for pep in fixed_peptides:
            self.passedPeptides.append(pep)

    def filter_by_name(self, df, output):
        df = pd.read_csv(df, sep='\t')
        names = self.transcriptomePeptides["proteinIds"].tolist()
        names = self.genomePeptides["proteinIds"].tolist()
        fixednames = []
        for name in names:
            splat = name.split(",")
            for n in splat:
                fixednames.append(n)
        df = df[df["Final Entries"].isin(fixednames)]
        df.to_csv(output, sep='\t', index=False)

    def filter_results(self, df, output):
        df = pd.read_csv(df, sep='\t')
        print(self.passedPeptides)
        df = df[df["Fixed Peptides"].isin(self.passedPeptides)]
        df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: <genome_pep_cat> <transcriptome_pep_cat> <df_results> <output_filtered_df>\n'
              'filters the results using a peptide FDR of 0.01')
    else:
        data = PeptideFDR(genome_pep_cat=sys.argv[1], transcriptome_pep_cat=sys.argv[2])
        # data.filter_results(df=sys.argv[3], output=sys.argv[4])
        data.filter_by_name(df=sys.argv[3], output=sys.argv[4])
