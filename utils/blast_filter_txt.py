import os
import sys

from Bio import SeqIO


class BlastFilterTXT:
    def __init__(self, results, fasta, evalue=0.01, identity=90, score=50, cov=90):
        self.results = results
        self.fasta = fasta
        self.eValue = evalue
        self.identity = identity
        self.score = score
        self.cov = cov

        self.toRemove = []

    def parse(self):
        with open(self.results, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if not line.startswith("#"):
                    cols = line.split("\t")
                    identity = float(cols[2])
                    score = float(cols[11])
                    evalue = float(cols[10])
                    query = cols[0]
                    cov = float(cols[12])
                    if identity >= self.identity:
                        if score >= self.score:
                            if cov >= self.cov:
                                if evalue <= self.eValue:
                                    self.toRemove.append(query)

    def filter(self, output):
        records = SeqIO.parse(self.fasta, 'fasta')
        filtered = []
        for record in records:
            entry = str(record.description)
            if entry not in self.toRemove:
                filtered.append(f'>{entry}\n{str(record.seq)}\n')
        with open(output, 'w') as handler:
            handler.writelines(filtered)

if __name__ == '__main__':
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("usage: .py <blast_results_outfmt_7_with_qcovs> <fasta_to_filter> <output>")
    else:
        data = BlastFilterTXT(results=sys.argv[1], fasta=sys.argv[2])
        data.parse()
        data.filter(output=sys.argv[3])