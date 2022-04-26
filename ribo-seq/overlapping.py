import sys
import os

import pandas as pd


class Overlapping(object):
    def __init__(self, slopped_intersected_gtf, intersected_gtf):
        self.intersected = intersected_gtf
        self.sloppedIntersected = slopped_intersected_gtf

        self.intersectedNovelORFs = {}

    def check_intersected(self, gtf, slopped=False):

        with open(gtf, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split("\t")

                # ref stuff
                ref = cols[:9]
                ref_attrs = ref[8].split(";")
                gene = ''
                transcript = ''
                gene_name = ''
                for a in ref_attrs:
                    if a.startswith(" "):
                        a = a[1:]
                    if 'gene_id' in a:
                        gene = self.__get_attribute_id(a)
                    if 'transcript_id' in a:
                        transcript = self.__get_attribute_id(a)
                    if 'gene_name' in a:
                        gene_name = self.__get_attribute_id(a)



                # novel stuff
                novel = cols[9:]
                novel_attrs = novel[8].split(";")
                novel_orf = ''
                for a in novel_attrs:
                    if 'orf_id' in a:
                        novel_orf = a.split(" ")[1].replace("\"", "")
                if novel_orf not in self.intersectedNovelORFs:
                    self.intersectedNovelORFs[novel_orf] = {'ref_gene_id': [], 'ref_transcript': [],
                                                            'ref_gene_name': [], 'surrounding_gene_ids': [],
                                                            'surrounding_transcripts': [], 'surrounding_gene_names': []}
                if not slopped:
                    self.intersectedNovelORFs[novel_orf]['ref_gene_id'].append(gene)
                    self.intersectedNovelORFs[novel_orf]['ref_transcript'].append(transcript)
                    self.intersectedNovelORFs[novel_orf]['ref_gene_name'].append(gene_name)
                elif slopped:

                    if gene not in self.intersectedNovelORFs[novel_orf]['ref_gene_id']: # we dont wanna say that a gene is surrounding an ORF if this same one
                        self.intersectedNovelORFs[novel_orf]['surrounding_gene_ids'].append(gene) # is already overlapping
                    else:
                        self.intersectedNovelORFs[novel_orf]['surrounding_gene_ids'].append('not found')

                    if transcript not in self.intersectedNovelORFs[novel_orf]['ref_transcript']:
                        self.intersectedNovelORFs[novel_orf]['surrounding_transcripts'].append(transcript)
                    else:
                        self.intersectedNovelORFs[novel_orf]['surrounding_transcripts'].append('not found')

                    if gene_name not in self.intersectedNovelORFs[novel_orf]['ref_gene_name']:
                        self.intersectedNovelORFs[novel_orf]['surrounding_gene_names'].append(gene_name)
                    else:
                        self.intersectedNovelORFs[novel_orf]['surrounding_gene_names'].append('not found')

                # print(ref)
                # print(novel)
                # break


    @staticmethod
    def __get_attribute_id(attribute):
        return attribute.split(" ")[1].replace("\"", "")


    def save(self, output):
        novel = []
        ref_gene_id = []
        ref_transcript = []
        ref_gene_name = []
        surr_gene_id = []
        surr_transcripts = []
        surr_gene_names = []

        for smorf in self.intersectedNovelORFs:
            for i in range(len(self.intersectedNovelORFs[smorf]['ref_gene_id'])):
                novel.append(smorf)
                ref_gene_id.append(''.join(self.intersectedNovelORFs[smorf]['ref_gene_id'][i]))
                ref_transcript.append(''.join(self.intersectedNovelORFs[smorf]['ref_transcript'][i]))
                ref_gene_name.append(''.join(self.intersectedNovelORFs[smorf]['ref_gene_name'][i]))
            # for i in range(len(self.intersectedNovelORFs[smorf]['surrounding_gene_ids'])):
                surr_gene_id.append(''.join(self.intersectedNovelORFs[smorf]['surrounding_gene_ids'][i]).replace("not found", ''))
                surr_transcripts.append(''.join(self.intersectedNovelORFs[smorf]['surrounding_transcripts'][i]).replace("not found", ''))
                surr_gene_names.append(''.join(self.intersectedNovelORFs[smorf]['surrounding_gene_names'][i]).replace("not found", ''))
        data = {'novel_ORF': novel, 'overlapping_gene_id': ref_gene_id, 'overlapping_transcript': ref_transcript,
                'overlapping_gene_name': ref_gene_name, 'surrounding_gene_ids': surr_gene_id,
                'surrounding_transcripts': surr_transcripts, 'surrounding_gene_names': surr_gene_names}
        df = pd.DataFrame(data=data)
        df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    # data = Overlapping(slopped_intersected_gtf='/media/eduardo/DATA/Eduardo/adenovirus/updated_annotation_files/102/custom_for_bedtools/overlapping_orfs_ribocode_102_slopped',
    #                    intersected_gtf='/media/eduardo/DATA/Eduardo/adenovirus/updated_annotation_files/102/custom_for_bedtools/overlapping_orfs_ribocode_102')
    # data.check_intersected(data.intersected)
    # data.check_intersected(data.sloppedIntersected, slopped=True)
    # data.save(output='/media/eduardo/DATA/Eduardo/adenovirus/updated_annotation_files/102/custom_for_bedtools/overlapping_orfs.xls')

    # folder = '/media/eduardo/DATA/Eduardo/adenovirus/updated_annotation_files/421/custom_for_bedtools'
    # data_421 = Overlapping(slopped_intersected_gtf=f'{folder}/intersected_421_orfs_slopped',
    #                        intersected_gtf=f'{folder}/intersected_421_orfs_not_slopped')
    # data_421.check_intersected(data_421.intersected)
    # data_421.check_intersected(data_421.sloppedIntersected, slopped=True)
    # data_421.save(output=f'{folder}/overlapping_orfs_421.xls')
    if sys.argv[1] == '-h':
        print('usage: overlapping.py <slopped_intersected_gtf> <non_slopped_intersected_gtf> <output>')
    else:
        data = Overlapping(slopped_intersected_gtf=sys.argv[1], intersected_gtf=sys.argv[2])
        data.check_intersected(data.intersected)
        data.check_intersected(data.sloppedIntersected, slopped=True)
        data.save(output=sys.argv[3])