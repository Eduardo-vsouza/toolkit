import os

import pandas as pd
from Bio import SeqIO
import argparse


def get_args():
    parser = argparse.ArgumentParser(description="Filter fasta file based on blastP results. Removes annotated"
                                                 " proteins.")
    parser.add_argument("--blastResults", help="Blastp tab-delimited results in outfmt 6 format.")
    parser.add_argument("--eValue", default=0.01, type=float)
    parser.add_argument("--identity", default=0.95, type=float)
    parser.add_argument("--bitScore", default=50, type=float)
    parser.add_argument("--fasta", help="fasta file to be filtered.")
    parser.add_argument("--outdir", help="Output file name.")
    parser.add_argument("--proteome", help="Reference proteome to search for annotated proteins.")
    parser.add_argument("--outfile")
    args = parser.parse_args()
    return args


class BlasterP(object):
    def __init__(self, new_smorfs, ref_proteomes, outdir):
        self.smorfs = new_smorfs
        self.referenceProteome = ref_proteomes
        self.outdir = outdir
        self.__check_dir()

    def __check_dir(self):
        if not os.path.exists(self.outdir):
            os.system(f'mkdir {self.outdir}')

    def blast(self):
        cmd = f'blastp -query {self.smorfs} -subject {self.referenceProteome} -outfmt 6 -out {self.outdir}/blasted_to_reference.txt'
        os.system(cmd)


class BlastParser(object):
    def __init__(self, blast_results, outdir, outfile):
        self.outdir = outdir
        self.df = pd.read_csv(blast_results, sep='\t')
        self.df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                           'send', 'evalue', 'bitscore']
        self.__check_dir()
        self.outfile = outfile
    def __check_dir(self):
        if not os.path.exists(self.outdir):
            os.system(f'mkdir {self.outdir}')

    def filter_evalue(self, evalue):
        self.df = self.df[self.df["evalue"] <= evalue]

    def filter_identity(self, identity):
        self.df = self.df[self.df["pident"] >= identity]

    def filter_bitscore(self, bitscore):
        self.df = self.df[self.df["bitscore"] >= bitscore]

    def filter_fasta(self, fasta, outdir):
        entries = self.df["qseqid"].tolist()
        records = SeqIO.parse(fasta, 'fasta')
        filtered = []
        checker = []
        for record in records:
            if str(record.id) not in entries and str(record.id) not in checker:
                filtered.append(f'>{str(record.id)}\n{str(record.seq)}\n')
                checker.append(str(record.id))
        with open(f'{outdir}/{self.outfile}', 'w') as out:
            out.writelines(filtered)


if __name__ == '__main__':
    args = get_args()

    # blast = BlasterP(new_smorfs=args.fasta, ref_proteomes=args.proteome, outdir=args.outdir)
    # blast.blast()

    data = BlastParser(blast_results=args.blastResults, outdir=args.outdir, outfile=args.outfile)
    data.filter_evalue(evalue=args.eValue)
    data.filter_identity(identity=args.identity)
    data.filter_bitscore(bitscore=args.bitScore)
    data.filter_fasta(fasta=args.fasta, outdir=args.outdir)

