import os
import sys

import argparse

class AssemblyQuantifier(object):
    def __init__(self, merged_gtf, illumina_reads_folder, nanopore_1, nanopore_2, outdir, threads):
        self.gtf = merged_gtf
        self.illuminaReadsFolder = illumina_reads_folder
        self.nanoporeAlignments = [nanopore_1, nanopore_2]
        self.outdir = outdir
        self.__check_dir(self.outdir)
        self.__check_dir(f'{self.outdir}/abundances')
        self.__check_dir(f'{self.outdir}/gtf')
        self.threads = threads

    @staticmethod
    def __check_dir(folder):
        if not os.path.exists(folder):
            os.system(f'mkdir {folder}')

    def quantify_illumina(self):
        files = os.listdir(self.illuminaReadsFolder)

        for file in files:
            if file.endswith(".bam"):
                print(f'Quantifying transcripts for {file}\n')
                full = f'{self.illuminaReadsFolder}/{file}'
                cmd = f'stringtie -B -p {self.threads} --rf -G {self.gtf} -A {self.outdir}/abundances/{file} -o {self.outdir}/gtf/{file} ' \
                      f'{full}'
                os.system(cmd)
                print(f'Done quantifying transcripts for {file}\n')

    def quantify_nanopore(self):
        i = 0
        for file in self.nanoporeAlignments:
            i += 1
            name = f'{file.split("/")[-1]}_{i}'
            print(f'Quantifying transcripts for {name}\n')
            cmd = f'stringtie -B -e -p {self.threads} -G {self.gtf} -A {self.outdir}/abundances/{name} -o {self.outdir}/gtf/{name} ' \
                  f'{file}'
            os.system(cmd)
            print(f'Done quantifying transcripts for {name}\n')

def get_args():
    parser = argparse.ArgumentParser(description='Creates GTF abundances files using the merged GTF file. ')
    parser.add_argument('--gtf', help="The merged GTF file generated using StringTie --merge")
    parser.add_argument('--illumina_reads', help="The folder containing illumina alignments (sorted bam files)")
    parser.add_argument('--nanopore1', help="The first nanopore bam alignments")
    parser.add_argument('--nanopore2', help="The second nanopore bam alignments")
    parser.add_argument('--outdir')
    parser.add_argument('--threads')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    data = AssemblyQuantifier(merged_gtf=args.gtf, illumina_reads_folder=args.illumina_reads, nanopore_1=args.nanopore1,
                              nanopore_2=args.nanopore2, outdir=args.outdir, threads=args.threads)
    data.quantify_illumina()
    data.quantify_nanopore()


