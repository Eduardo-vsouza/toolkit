import sys

import pandas as pd
from Bio import SeqIO
from difflib import SequenceMatcher


class TranscriptRenamer(object):
    def __init__(self, fasta_x, fasta_y):
        self.fastaX = fasta_x
        self.fastaY = fasta_y

        self.transcriptsX, self.entriesX = self.__extract_sequences(self.fastaX)
        self.transcriptsY, self.entriesY = self.__extract_sequences(self.fastaY)

        self.outputReplaced = f'{self.fastaX.replace(".fasta", "")}_replaced.fasta'


    @staticmethod
    def __extract_sequences(fasta):
        transcripts = {}
        entries = []
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            desc = str(record.description)
            if desc not in transcripts:
                entries.append(desc)
                transcripts[str(record.seq)] = desc
        return transcripts, entries

    def compare(self):
        for rna in self.transcriptsX:
            if 'MSTRG' in self.transcriptsX[rna]:
                self.transcriptsX[rna] = f'{self.transcriptsX[rna]}_alt'
        new_x = []
        replaced = {}
        to_change = []  # if a new renamed entry was already in the fasta file, change it to another thing
        for rna in self.transcriptsX:
            if self.transcriptsX[rna] in to_change:
                name = f'{self.transcriptsX[rna]}_alt'
            else:
                name = self.transcriptsX[rna]
            if 'MSTRG' in name:
                if rna in self.transcriptsY:
                    # if SequenceMatcher(None, name, self.transcriptsY[rna]).ratio() >= 0.9:
                    if name == self.transcriptsY[rna]:
                        # to_change.append(name)
                        new_x.append(f'>{name}\n{rna}\n')
                    else:
                        new_x.append(f'>{self.transcriptsY[rna]}\n{rna}\n')
                        replaced[name] = self.transcriptsY[rna]
                        # if self.transcriptsY[rna] in self.entriesX:
                        #     to_change.append(self.transcriptsY[rna])

                else:
                    new_x.append(f'>{name}\n{rna}\n')
                    # to_change.append(name)
            else:
                new_x.append(f'>{name}\n{rna}\n')
                # to_change.append(name)
                # replaced[self.transcriptsX[rna]] = self.transcriptsX[rna]

        data = {'original': list(replaced.keys()), 'replaced': list(replaced.values())}
        df = pd.DataFrame(data)
        df.to_csv(f'{self.fastaX.replace(".fasta", "")}_changed_names.txt', sep='\t', index=False)
        self.outputReplaced = f'{self.fastaX.replace(".fasta", "")}_replaced.fasta'
        with open(f'{self.fastaX.replace(".fasta", "")}_replaced.fasta', 'w') as out:
            out.writelines(new_x)

    def check_duplicates(self):
        duplicated = []
        checker = []
        records = SeqIO.parse(self.outputReplaced, 'fasta')
        for record in records:
            name = str(record.description)
            if name not in checker:
                checker.append(name)
            else:
                duplicated.append(name)
        print('These entries are duplicated: ', duplicated)

    def check_sequence_duplicates(self):
        duplicated = []
        checker = []
        records = SeqIO.parse(self.outputReplaced, 'fasta')
        for record in records:
            seq = str(record.seq)
            if seq not in checker:
                checker.append(seq)
            else:
                duplicated.append(seq)
        print(duplicated)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: transcript_renamer.py <fasta_x> <fasta_y>\nfasta_x is the one that is going to have their '
              'transcript names replaced by the ones in fasta_y')
    else:
        data = TranscriptRenamer(fasta_x=sys.argv[1], fasta_y=sys.argv[2])
        data.compare()
        # data.check_duplicates()
        # data.check_sequence_duplicates()
