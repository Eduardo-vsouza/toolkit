import sys

import pandas as pd
from Bio import SeqIO


class CDSExtractor(object):
    def __init__(self, results, genome):
        self.results = pd.read_csv(results, sep='\t')
        self.genome = self.__get_genome_sequence(genome)


    @staticmethod
    def __get_genome_sequence(genome):
        seq = []
        records = SeqIO.parse(genome, 'fasta')
        for record in records:
            seq.append(str(record.seq))
        return seq[0]

    def extract_nucleotides(self, output, subset='genome'):
        sequences = []
        starts = []
        coords = self.results["Genome Coordinates"].tolist()
        entries = self.results["Final Entries"].tolist()
        for entry, coord in zip(entries, coords):
            splat = coord.split("-")
            if coord != 'not found':
                if subset == 'genome':
                    if 'reverse' in entry:
                        start, end = int(splat[1]), int(splat[0])
                        cds = self.__get_reverse_complement(start, end - 3)
                    else:
                        start, end = int(splat[0]), int(splat[1])
                        cds = self.__get_forward_cds(start, end)
                elif subset == 'transcriptome':
                    if int(splat[0]) > int(splat[1]):
                        start, end = int(splat[1]), int(splat[0])
                        cds = self.__get_reverse_complement(start, end)
                    else:
                        start, end = int(splat[0]), int(splat[1])
                        cds = self.__get_forward_cds(start, end)
                else:
                    break
                sequences.append(cds)
                starts.append(cds[:3])
            else:
                sequences.append('not found')
                starts.append('not found')
        self.results.insert(10, "orf_sequence_nucleotides", sequences)
        self.results.insert(11, "start_codon", starts)
        self.results.to_csv(output, sep='\t', index=False)

    def __get_forward_cds(self, start, end):
        return self.genome[start-1: end]

    def __get_reverse_complement(self, start, end):
        """ start and end must disregard strand position. """
        rev = self.genome[start-1: end][::-1]
        nucs = {'A': 'T', 'T': 'A',
                'G': 'C', 'C': 'G'}
        comp = ''.join([nucs[i] for i in rev])
        return comp

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: <genome> <reslts> <output> <subset>\n'
              'genome: genome fasta file\n'
              'subset: genome or transcriptome\n')
    else:
        data = CDSExtractor(genome=sys.argv[1], results=sys.argv[2])
        data.extract_nucleotides(output=sys.argv[3], subset=sys.argv[4])

