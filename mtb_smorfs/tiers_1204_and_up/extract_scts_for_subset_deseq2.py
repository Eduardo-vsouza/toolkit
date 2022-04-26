import os
import sys

from Bio import SeqIO


class Extractor:
    def __init__(self, fasta):
        self.fasta = fasta

    def extract(self, output):
        names = []
        records = SeqIO.parse(self.fasta, 'fasta')
        for record in records:
            entry = str(record.description)
            if 'tORF' in entry:
                names.append(f'{entry}\n')

        with open(output, 'w') as handler:
            handler.writelines(names)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: <fasta> <output>')
    else:
        data = Extractor(fasta=sys.argv[1])
        data.extract(output=sys.argv[2])
