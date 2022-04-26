import sys

from Bio import SeqIO
import pandas as pd


class SCTRenamer:
    def __init__(self, transcripts_fasta, scts_table):
        """

        :param transcripts_fasta: transcript file generated with gffread after performing the assembly with Stringtie
        :param scts_table: table generated with scts_extractor.py
        """
        self.transcripts = transcripts_fasta
        self.df = pd.read_csv(scts_table, sep='\t')

        self.SCTs = {}
        self.__get_scts()

    def __get_scts(self):
        scts = self.df["scts"].tolist()
        torfs = self.df["torfs"].tolist()
        for sct, torf in zip(scts, torfs):
            self.SCTs[sct] = torf

    def rename(self, output):
        records = SeqIO.parse(self.transcripts, 'fasta')
        fasta = []
        for record in records:
            entry = str(record.description).replace("gene-", "")
            if entry in self.SCTs:
                new = self.SCTs[entry]
                fasta.append(f'>{new}\n{str(record.seq)}\n')
            else:
                fasta.append(f'>{entry}\n{str(record.seq)}\n')
        with open(output, 'w') as handler:
            handler.writelines(fasta)


if __name__ == '__main__':
    data = SCTRenamer(transcripts_fasta=sys.argv[1], scts_table=sys.argv[2])
    data.rename(output=sys.argv[3])