import os
import sys

from Bio import SeqIO



def remove_dupli(folder, outdir):
    files = os.listdir(folder)
    for file in files:
        if file.endswith('fasta'):
            checker = []
            new = []
            records = SeqIO.parse(f'{folder}/{file}', 'fasta')
            for record in records:
                seq = str(record.seq)
                if seq not in checker:
                    checker.append(seq)
                    new.append(f'>{str(record.description)}\n{seq}\n')
            with open(f'{outdir}/{file.replace(".fasta", "")}_deduplicated.fasta', 'w') as handler:
                handler.writelines(new)


if __name__ == '__main__':
    remove_dupli(folder=sys.argv[1], outdir=sys.argv[2])