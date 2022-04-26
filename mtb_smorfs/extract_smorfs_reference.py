import sys

from Bio import SeqIO


class Smorfetcher(object):
    def __init__(self, fasta):
        self.fasta = fasta

    def get_smorfs(self, output):
        records = SeqIO.parse(self.fasta, 'fasta')
        fasta = []
        checker = []
        for record in records:
            seq = str(record.seq)
            if seq not in checker:
                if len(seq) <= 100:
                    fasta.append(record)
        SeqIO.write(fasta, output, 'fasta')


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: extract_smorfs_reference.py <input> <output>')
    else:
        data = Smorfetcher(fasta=sys.argv[1])
        data.get_smorfs(output=sys.argv[2])