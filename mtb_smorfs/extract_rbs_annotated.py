#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO


class UpstreamAnalysis(object):
    def __init__(self, gtf, genome):
        self.gtf = gtf
        self.genome = self.__extract_genome_sequence(genome)

        self.coordinates = {}
        self.cds = {}
        self.startCodons = {}
        self.upstream = {}

    @staticmethod
    def __extract_genome_sequence(genome):
        seq = []
        records = SeqIO.parse(genome, 'fasta')
        for record in records:
            seq.append(str(record.seq))
        return seq[0]

    def extract_coordinates(self):
        with open(self.gtf, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                cols = line.split("\t")
                feature = cols[2]
                if feature == 'CDS':
                    start, end = int(cols[3]), int(cols[4])
                    if end - start <= 300:
                        attrs = cols[8].split(";")
                        strand = cols[6]
                        for a in attrs:
                            if 'gene_id' in a:
                                gene = a.split(" ")[2].replace("\"", "")
                                self.coordinates[gene] = (start, end, strand)

    def extract_cds(self):
        for gene in self.coordinates:
            strand = self.coordinates[gene][2]
            if strand == '+':
                start = self.coordinates[gene][0]
                end = self.coordinates[gene][1]
                cds = self.genome[start-1: end-3]
            else:
                start = self.coordinates[gene][0]
                end = self.coordinates[gene][1]
                cds = self.__get_reverse_complement(start, end)

            self.cds[gene] = cds
            self.startCodons[gene] = cds[:3]

    def extract_upstream(self):
        for gene in self.coordinates:
            strand = self.coordinates[gene][2]
            start = self.coordinates[gene][0]
            end = self.coordinates[gene][1]
            if strand == '+':
                upstream = self.genome[start-23: start-1]
            else:
                upstream = self.__get_reverse_complement(end+1, end+22)
            self.upstream[gene] = upstream

    def __get_reverse_complement(self, start, end):
        """ start and end must disregard strand position. """
        rev = self.genome[start-1: end][::-1]
        nucs = {'A': 'T', 'T': 'A',
                'G': 'C', 'C': 'G'}
        comp = ''.join([nucs[i] for i in rev])
        return comp

    def save(self, output):
        genes = []
        upstream = []
        starts = []
        cds = []
        for gene in self.upstream:
            genes.append(gene)
            upstream.append(self.upstream[gene])
            starts.append(self.startCodons[gene])
            cds.append(self.cds[gene])
        df = pd.DataFrame(data={'gene': genes, 'upstream_seq_22': upstream, 'cds': cds, 'start_codon': starts})
        df.to_csv(output, sep='\t', index=False)




if __name__ == '__main__':
    data = UpstreamAnalysis(genome='/media/eduardo/New Volume/Eduardo/mtb_annotation_files/genome.fasta',
                            gtf='/media/eduardo/New Volume/Eduardo/mtb_annotation_files/mtb.gtf')
    data.extract_coordinates()
    data.extract_cds()
    data.extract_upstream()
    data.save('/media/eduardo/New Volume/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info.txt')
    # print(len(data.cds))
    # print(data.cds['gene-Rv0033'])
    # print(data.startCodons)