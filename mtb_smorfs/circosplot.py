import os
import sys
import collections

import pycircos
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


Garc = pycircos.Garc
Gcircle = pycircos.Gcircle


class Circos:
    def __init__(self, df, tiers=('T1', 'T2', 'T3')):
        self.circle = Gcircle()
        self.tiers = tiers
        self.df = pd.read_csv(df, sep='\t')
        self.df = self.df[self.df["Genome Coordinates"] != "Genome Coordinates"]
        self.df = self.df.drop_duplicates(subset=["Final Entries"])
        self.df = self.df[self.df["Tier"].isin(tiers)]
        self.coordinates = collections.defaultdict(dict)
        self.orfDict = {}

        self.outerMostCircleRanges = (940, 990) # original 950, 1000
        self.outerCircleTwoRanges = (870, 920)

    def create_circle(self):
        print('Creating circles')
        name = 'chr1'
        length = 4411532
        arc = Garc(arc_id=name, size=length, interspace=1, raxis_range=self.outerMostCircleRanges, labelposition=60, label_visible=False,
                   facecolor='#ffffff')
        self.circle.add_garc(arc)
        self.circle.set_garcs()

        # arc = Garc(arc_id='chr2', size=length, interspace=1, raxis_range=self.outerCircleTwoRanges, labelposition=60, label_visible=False,
        #            facecolor='#ffffff', linewidth=0.5)
        # self.circle.add_garc(arc)
        # self.circle.set_garcs()

        # self.circle.set_garcs(-65, 245)

    def add_orfs(self):
        print('Adding novel ORFs')
        coordinates = self.df["Genome Coordinates"].tolist()
        orfs = self.df["Final Entries"].tolist()

        for coord, name in zip(coordinates, orfs):
            orf = name
            name = 'chr1'
            splat = coord.split("-")
            start, end = int(splat[0]), int(splat[1])
            if start > end:
                start, end = end, start
            if name not in self.coordinates:
                self.coordinates[name] = {}
                self.coordinates[name]['start'] = []
                self.coordinates[name]['width'] = []
                self.coordinates[name]['color'] = []

            self.coordinates[name]['start'].append(start)
            self.coordinates[name]['width'].append((end-start)+500)

            # defines different colors for gORFs and tORFs
            if 'tORF' in orf:
                color = '#00af2a'
            else:
                color = "#0083ff"
            # color = '#e81d1d'
            self.coordinates[name]['color'].append(color)
            entry_splat = orf.split("_")
            orf_name = f'{entry_splat[0]}_{entry_splat[2]}'
            self.orfDict[orf_name] = {'start': start, 'width': (end-start)+500}

        # adds orfs from dictionary
        for key in self.coordinates:
            # print(self.coordinates[key]['start'])
            # print(self.circle._garc_dict['chr1'])
            self.circle.barplot(key, data=[1]*len(self.coordinates[key]["start"]),
                                positions=self.coordinates[key]["start"],
                                width=self.coordinates[key]["width"], raxis_range=self.outerMostCircleRanges,
                                facecolor=self.coordinates[key]["color"])

    def add_annotated(self, anno_color):
        print('Adding annotated')
        annotated = {}
        for transcript in self.transcriptData:
            # for coord, name in zip(coordinates, orfs):
            name = 'chr1'
            start = self.transcriptData[transcript]['start']
            end = self.transcriptData[transcript]['end']
            if end - start <= 300:
                color = '#3a70ce'
            else:
                color = anno_color
            # start, end = int(splat[0]), int(splat[1])
            if start > end:
                start, end = end, start
            if name not in annotated:
                annotated[name] = {}
                annotated[name]['start'] = []
                annotated[name]['width'] = []
                annotated[name]['color'] = []

            annotated[name]['start'].append(start)
            annotated[name]['width'].append((end - start) + 500)

            # defines different colors for gORFs and tORFs
            # if 'tORF' in orf:
            #     color = '#00af2a'
            # else:
            #     color = "#0083ff"
            # color = '#e81d1d'
            annotated[name]['color'].append(color)

            # adds orfs from dictionary
        for key in annotated:
            # print(self.coordinates[key]['start'])
            # print(self.circle._garc_dict['chr1'])
            self.circle.barplot(key, data=[1] * len(annotated[key]["start"]),
                                positions=annotated[key]["start"],
                                width=annotated[key]["width"], raxis_range=self.outerCircleTwoRanges,
                                facecolor=annotated[key]["color"])

    def add_co_expression_links(self, corr_table, gtf_annotation_table):
        links = {}
        print('Adding co-expression links')
        co_expressed = self.__get_co_expressed(corr_table)
        df = self.df.drop_duplicates(subset=['Final Entries'])
        entries = df["Final Entries"].tolist()
        coordinates = df["Genome Coordinates"].tolist()

        transcript_data = self.__get_transcripts_from_gtf_table(gtf_annotation_table)

        for gene in co_expressed:
            for i, entry in enumerate(entries):
                if 'tORF' in entry:
                    splat = entry.split("_")
                    torf = f'{splat[0]}_{splat[2]}'   # gets the tORF number, like tORF_20324, from the entry
                    if torf == gene:                # important because this is how they are in the co-expression nets
                        coords = coordinates[i].split("-")
                        start, end = int(coords[0]), int(coords[1])
                        if start > end:
                            start, end = end, start
                        source = ('chr1', start, end, 630)
                        for target in co_expressed[gene]:  # we need the same as 'source' for the co expressed genes
                            if 'tORF' not in target:
                                target_start, target_end = transcript_data[target]['start'], transcript_data[target]['end']
                                destination = ('chr1', target_start, target_end, 630)
                                self.circle.chord_plot(source, destination, facecolor='#89ffb5')

    def __get_transcripts_from_gtf_table(self, gtf_table):
        transcript_data = {}
        df = pd.read_csv(gtf_table, sep='\t')
        transcripts = df["transcript_id"].tolist()
        starts = df["start"].tolist()
        ends = df["end"].tolist()
        for i, transcript in enumerate(transcripts):
            name = transcript.replace("gene-", "")  # the co-expression table doesn't have the 'gene-' prefix for genes
            transcript_data[name] = {'start': int(starts[i]), 'end': int(ends[i])}
        self.transcriptData = transcript_data
        return transcript_data


    def __get_co_expressed(self, corr_table):
        co_expressed = {}
        df = pd.read_csv(corr_table, sep='\t')
        df = df[df["gene1"].str.contains('tORF')]
        genes = df["gene1"].tolist()
        targets = df["gene2"].tolist()
        for gene, target in zip(genes, targets):
            if gene not in co_expressed:
                co_expressed[gene] = []
            co_expressed[gene].append(target)
        return co_expressed

    def add_transcript_expression(self, counts_table):
        print('Adding transcript expression')
        tpms = {}
        arcdata_dict = collections.defaultdict(dict)
        values_all = []
        df = pd.read_csv(counts_table, sep='\t')
        genes = df["gene"].tolist()
        control1 = df["control_1"].tolist()
        control2 = df["control_2"].tolist()
        control3 = df["control_3"].tolist()
        # print(self.transcriptData)
        for i, gene in enumerate(genes):
            gene = gene.replace("gene-", "")
            mean = np.log1p(np.mean([control1[i], control2[i], control3[i]]))
            # print(mean)
            if gene not in tpms:
                tpms[gene] = {}
            tpms[gene]['tpm'] = mean
            tpms[gene]['start'] = self.transcriptData[gene]['start']
            tpms[gene]['end'] = self.transcriptData[gene]['end']
        arcdata_dict['chr1']["positions"] = []
        arcdata_dict['chr1']["widths"] = []
        arcdata_dict['chr1']["values"] = []
        for gene in tpms:
            # if gene not in arcdata_dict:

            # print(gene, tpms[gene]['start'], tpms[gene]['tpm'])
            arcdata_dict['chr1']["positions"].append(tpms[gene]['start'])
            arcdata_dict['chr1']["widths"].append(tpms[gene]['end']-tpms[gene]['start'])
            arcdata_dict['chr1']["values"].append(tpms[gene]['tpm'])
            values_all.append(tpms[gene]['tpm'])

        vmin, vmax = min(values_all), max(values_all)
        # print(vmin, vmax)
        # for key in arcdata_dict:
        #     self.circle.heatmap(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
        #                    width=arcdata_dict[key]["widths"], raxis_range=[700, 870], vmin=vmin, vmax=vmax,
        #                    cmap=plt.cm.viridis)
        for key in arcdata_dict:
            # print(arcdata_dict[key]["values"])
            self.circle.barplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
                           width=arcdata_dict[key]["widths"], base_value=0.0,
                           rlim=[vmin-0.05*abs(vmin), vmax+0.05*abs(vmax)],
                           raxis_range=[680, 800], facecolor="#c00000", spine=False)

    def add_conservation(self, color):
        blastp_myco = '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/phylogenetics/blastp/mycobacteria/conservation_by_protein'
        blastp_actino = '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/phylogenetics/blastp/actinobacteria/conservation_by_protein'
        blastp = [blastp_myco, blastp_actino]
        tblastn_myco = '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/phylogenetics/tblastn/mycobacteria/conservation_by_protein'
        tblastn_actino = '/media/eduardo/gold/Eduardo/smorfs_mtb_090222/phylogenetics/tblastn/actinobacteria/conservation_by_protein'
        tblastn = [tblastn_myco, tblastn_actino]

        def get_conservation(blast_list):
            conservation = {}
            for folder in blast_list:
                files = os.listdir(folder)
                for file in files:
                    tier = file.split("_")[1]
                    if tier in self.tiers:
                        df = pd.read_csv(f'{folder}/{file}', sep='\t')
                        proteins = df["protein"].tolist()
                        ortologous = df["ortologous"].tolist()
                        for protein, orto in zip(proteins, ortologous):
                            splat = protein.split("_")
                            name = f'{splat[0]}_{splat[2]}'
                            if name not in conservation:
                                conservation[name] = 0
                            conservation[name] += orto
            return conservation

        blastp_conservation = get_conservation(blastp)
        tblastn_conservation = get_conservation(tblastn)
        arcdata_dict = collections.defaultdict(dict)
        values_all = []

        arcdata_dict['chr1']["positions"] = []
        arcdata_dict['chr1']["widths"] = []
        arcdata_dict['chr1']["values"] = []

        def add_to_circle(conservation):
            for gene in conservation:
                # if gene not in arcdata_dict:

                # print(gene, tpms[gene]['start'], tpms[gene]['tpm'])
                # print(gene)
                # print(self.orfDict[gene])
                arcdata_dict['chr1']["positions"].append(self.orfDict[gene]['start'])
                arcdata_dict['chr1']["widths"].append(self.orfDict[gene]['width'])
                arcdata_dict['chr1']["values"].append(blastp_conservation[gene])
                values_all.append(blastp_conservation[gene])

            vmin, vmax = min(values_all), max(values_all)
            for key in arcdata_dict:
                # self.circle.heatmap(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"],
                #                width=arcdata_dict[key]["widths"], raxis_range=[780, 870], vmin=vmin, vmax=vmax,
                #                cmap=plt.cm.RdBu)
                self.circle.barplot(key, data=np.log1p(arcdata_dict[key]["values"]), positions=arcdata_dict[key]["positions"],
                                    width=arcdata_dict[key]["widths"], base_value=0.0,
                                    rlim=np.log1p([vmin - 0.05 * abs(vmin), vmax + 0.05 * abs(vmax)]),
                                    raxis_range=[810, 860], facecolor=color, spine=True)
        add_to_circle(blastp_conservation)
        add_to_circle(tblastn_conservation)


if __name__ == '__main__':
    data = Circos(df='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/tier_analyses_1204/cat_transcriptome_genome_pep_fdr_tiers.csv')
    data.create_circle()
    data.add_co_expression_links(corr_table=f'/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/reads/isoniazid_GSE165581/quantification/all_together/Cor_table_filter.txt',
                                 gtf_annotation_table='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/assembly_chr_gene_transcript_start_end_gtf.xls')
    data.add_transcript_expression(counts_table='/media/eduardo/gold/Eduardo/smorfs_mtb_090222/transcriptomics_analysis/quantification/mtab-1616_for_deseqcounts.txt')
    data.add_annotated(anno_color='#a969d6')
    data.add_orfs()
    data.add_conservation('#ff5500')
    plt.show()
