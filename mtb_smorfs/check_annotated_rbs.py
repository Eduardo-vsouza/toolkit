import os
import subprocess

import pandas as pd
from Bio import SeqIO


class ShineDalgarno(object):
    def __init__(self, df, rrna):
        self.df = pd.read_csv(df, sep='\t')
        self.rRNA = self.__get_sequences(rrna)[::-1][:13]
        self.freeAlignPath = '/home/eduardo/programs/uproteins_1.0.3/uProteInS/dependencies/free2bind/free_align.pl'
        self.upstream = self.df["upstream_seq_22"].tolist()

    def add_rbs(self, output):
        rbs = []
        energies = []
        for upstream in self.upstream:
            cmd = f'{self.freeAlignPath} -e {upstream} {self.rRNA}'
            energy = subprocess.check_output(cmd, shell=True).strip().rstrip()
            freeEnergy = float(energy)
            sd_seq = self.__check_rbs(energy)
            energies.append(freeEnergy)
            rbs.append(sd_seq)
        self.df.insert(3, "shine_dalgarno", rbs)
        self.df.insert(4, "free_energy", energies)
        self.df.to_csv(output, sep='\t', index=False)

    @staticmethod
    def __get_sequences(fasta):
        records = SeqIO.parse(fasta, 'fasta')
        seqs = [str(record.seq) for record in records][0]
        return seqs

    @staticmethod
    def __check_rbs(rbs):
        if float(rbs) >= -3.4535:
            sd_seq = "SD absent"
        elif -8.4 < float(rbs) < -3.4535:
            sd_seq = "Present. Low to moderate binding."
        elif float(rbs) <= -8.4:
            sd_seq = "Present. Strong binding."
        else:
            sd_seq = "SD absent"
        return sd_seq


if __name__ == '__main__':
    folder = '/media/eduardo/New Volume/Eduardo/mtb_annotation_files/'
    data = ShineDalgarno(df='/media/eduardo/New Volume/Eduardo/mtb_annotation_files/gtf_smorfs_genome_nuc_info.txt',
                         rrna='/media/eduardo/New Volume/Eduardo/mtb_annotation_files/16s_rrna.fna')
    data.add_rbs(output=f'{folder}/gtf_smorfs_genome_nuc_info_with_rbs.txt')


