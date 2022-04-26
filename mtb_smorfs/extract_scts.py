import pandas as pd


class SCTs:
    def __init__(self, df, tier_list):
        self.df = pd.read_csv(df, sep='\t')
        self.df = self.df[self.df["Tier"].isin(tier_list)]

    def get_scts(self, output):
        scts = []
        full = []
        torfs = []
        entries = self.df["Final Entries"].tolist()
        for entry in entries:
            splat = entry.split("_")
            sct = splat[1]
            torf = f'{splat[0]}_{splat[2]}'
            scts.append(sct)
            full.append(entry)
            torfs.append(torf)
        df = pd.DataFrame(data={'scts': scts, 'full_entries': full, 'torfs': torfs})
        df = df.drop_duplicates(subset=['scts'])
        df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    # data = SCTs('/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/transcriptome_results_reformatted_0402_with_tiers.xls',
    #             tier_list=['T1', 'T2', 'T3'])
    # data.get_scts(output='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/torfs_scts_no_dupli.txt')

    "13/04/22"
    data = SCTs('/media/eduardo/gold/Eduardo/smorfs_mtb_090222/Transcriptome/Results/reanalysis_pep_fdr/transcriptome_results_with_tiers_pep_fdr_tiers_redone.xls',
                tier_list=['T1', 'T2', 'T3'])
    data.get_scts(output='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/torfs_scts_no_dupli.txt')
