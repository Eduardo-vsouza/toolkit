import pandas as pd


class Observer:
    def __init__(self, df, subset):
        self.subset = subset
        self.df = df

        self.subsetInfo = {}
        self.__get_info()

    def __get_info(self):
        if self.subset == 'ta':  # individual TA sites
            essential = self.df["Essentiality State"].tolist()
            coordinates = self.df["Coordinate"].tolist()
            self.subsetInfo['essentiality'] = essential
            self.subsetInfo['coordinates'] = coordinates

        elif self.subset == 'unannotated':
            start = self.df["Start Coordinate"].tolist()
            end = self.df["End Coordinate"].tolist()
            essential = self.df["Final Call"].tolist()
            infos = {'start': start, 'end': end, 'essential': essential}
            for i in infos:
                self.subsetInfo[i] = infos[i]

        elif self.subset == 'genomic_features':
            start = self.df["Start Coordinate"].tolist()
            end = self.df["End Coordinate"].tolist()
            essential = self.df["Final Call"].tolist()
            feature = self.df["Genomic Feature Type"].tolist()
            infos = {'start': start, 'end': end, 'essential': essential, 'feature': feature}
            for i in infos:
                self.subsetInfo[i] = infos[i]

class Essentiality(Observer):
    def __init__(self, df, subset):
        self.df = pd.read_csv(df, sep='\t')
        super().__init__(self.df, subset)

    def compare(self, cat_results_genome_transcriptome, tier_list, output):
        results = pd.read_csv(cat_results_genome_transcriptome, sep='\t')
        results = results[results["Final Entries"] != 'Final Entries']
        if self.subset == 'ta':
            smorfs_calls, smorfs_tiers = self.__compare_ta(results)
            smorfs = []
            calls = []
            tiers = []
            for i in smorfs_calls:
                if smorfs_calls[i][0] != 'none':
                    if smorfs_tiers[i] in tier_list:
                        print(set(smorfs_calls[i]))
                        print(i)
                        print(smorfs_tiers[i])
                        print('\n')
                        smorfs.append(i)
                        calls.append(f'{set(smorfs_calls[i])}'[1:-1].replace("\"", "")[1:-1])
                        tiers.append(smorfs_tiers[i])
            df = pd.DataFrame(data={'smorfs': smorfs, 'tier': tiers, 'essentiality': calls})
            df.to_csv(output, sep='\t', index=False)


        elif self.subset == 'genomic_features':
            smorfs_calls, smorfs_tiers, smorfs_features = self.__compare_features(results)
            print(smorfs_calls)
            for i in smorfs_calls:
                if smorfs_calls[i][0] != 'none':
                    if smorfs_tiers[i] in tier_list:
                        print(set(smorfs_calls[i]))
                        print(i)
                        print(smorfs_tiers[i])
                        print(smorfs_features[i])
                        print('\n')

    def __compare_features(self, results):
        starts = self.subsetInfo['start']
        ends = self.subsetInfo['end']
        call = self.subsetInfo['essential']
        features = self.subsetInfo['feature']

        smorfs_coords = results["Genome Coordinates"].tolist()
        entries = results["Final Entries"].tolist()
        tiers = results["Tier"].tolist()

        smorfs_calls = {}
        smorfs_tiers = {}
        smorfs_features = {}

        for i in range(len(smorfs_coords)):
            smorfs_tiers[entries[i]] = tiers[i]
            if entries[i] not in smorfs_calls:
                smorfs_calls[entries[i]] = []
                smorfs_features[entries[i]] = []
            for j in range(len(starts)):
                splat = smorfs_coords[i].split("-")  # smorfs
                start, end = starts[j], ends[j]  # sassetti coordinates
                smorf_start, smorf_end = int(splat[0]), int(splat[1])  # smorfs coordinates
                if smorf_start > smorf_end:
                    smorf_start, smorf_end = smorf_end, smorf_start
                if start in range(smorf_start, smorf_end) or end in range(smorf_start, smorf_end) or smorf_start in range(start, end) or smorf_end in range(start, end):
                    smorfs_calls[entries[i]].append(call[j])
                    smorfs_features[entries[i]].append(features[j])
            if len(smorfs_calls[entries[i]]) == 0:
                smorfs_calls[entries[i]].append('none')
        return smorfs_calls, smorfs_tiers, smorfs_features



    def __compare_ta(self, results):
        coords = self.subsetInfo['coordinates']
        call = self.subsetInfo['essentiality']
        smorfs_coords = results["Genome Coordinates"].tolist()
        entries = results["Final Entries"].tolist()
        tiers = results["Tier"].tolist()
        smorfs_calls = {}
        smorfs_tiers = {}
        i = 0
        for smorf_coord, smorf_entry in zip(smorfs_coords, entries):
            smorfs_tiers[smorf_entry] = tiers[i]
            i += 1
            if smorf_entry not in smorfs_calls:
                smorfs_calls[smorf_entry] = []
            for coord, call in zip(coords, call):
                splat = smorf_coord.split("-")
                start, end = int(splat[0]), int(splat[1])
                if start > end:
                    start, end = end, start
                if int(coord) in range(start, end):
                    smorfs_calls[smorf_entry].append(call)
            if len(smorfs_calls[smorf_entry]) == 0:
                smorfs_calls[smorf_entry].append('none')
        return smorfs_calls, smorfs_tiers




if __name__ == '__main__':
    folder = '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/sassetti_essentiality/cleaned/'
    data = Essentiality(df=f'{folder}/individual_ta_sites_essentiality.csv', subset='ta')
    # data.compare('/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls',
    #              tier_list=['T1', 'T2', 'T3'], output=f'/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/essentiality_tiers_1_to_3.xls')

    # features = Essentiality(df=f'{folder}/essentiality_calls_genomic_features.csv', subset='genomic_features')
    # features.compare('/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls',
    #              tier_list=['T1', 'T2', 'T3'])

    """ 14/03/22 """
    data.compare('/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tiers_analyses_0104/genome_transcriptome_cat_0402_with_tiers.xls',
                 tier_list=['T1', 'T2', 'T3'], output=f'/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/essentiality_tiers_1_to_3.xls')
