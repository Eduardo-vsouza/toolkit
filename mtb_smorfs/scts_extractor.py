import sys

import pandas as pd


class SCTHandling:
    def __init__(self, df, tier_list):
        self.df = pd.read_csv(df, sep='\t')
        self.tierList = tier_list
        self.__filter_tiers()

        self.SCTs = {}

    def __filter_tiers(self):
        self.df = self.df[self.df["Tier"].isin(self.tierList)]

    def get_scts(self):
        entries = self.df["Final Entries"].tolist()
        for entry in entries:
            splat = entry.split("_")
            transcript = splat[1].replace("gene-", "")
            if transcript not in self.SCTs:
                self.SCTs[transcript] = []
            orf = f'{splat[2]}'
            if orf not in self.SCTs[transcript]:
                self.SCTs[transcript].append(orf)

    def torfs_per_scts(self):
        for transcript in self.SCTs:
            if len(self.SCTs[transcript]) > 1:
                print(self.SCTs[transcript])

    def save_table(self, output):
        scts = []
        torfs = []
        for sct in self.SCTs:
            if len(self.SCTs[sct]) > 1:
                pattern = 'tORFs'
            else:
                pattern = 'tORF'
            torf_list = f'{pattern}_{"_".join(self.SCTs[sct])}'
            scts.append(sct)
            torfs.append(torf_list)
        data = {'scts': scts, 'torfs': torfs}
        df = pd.DataFrame(data=data)
        df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: .py <results_df> <comma_separated_tiers> <output_table>')
    else:
        data = SCTHandling(df=sys.argv[1], tier_list=sys.argv[2].split(","))
        data.get_scts()
        # print(data.SCTs)
        data.torfs_per_scts()
        data.save_table(output=sys.argv[3])