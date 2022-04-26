#!/usr/bin/python3

import os
import sys


def quantify_transcripts(reads_folder, salmon_index, output):
    """ Aligns short reads illumina data to a transcriptome fasta file (not to the genome) and quantifies it using
    salmon. """
    reads = os.listdir(reads_folder)
    for read in reads:
        if '_1.' in read and '.fastq' in read:
            read_1 = f'{reads_folder}/{read}'
            read_2 = f'{reads_folder}/{read.replace("_1.", "_2.")}'
            cmd = f'salmon quant -i {salmon_index} -l A -1 {read_1} -2 {read_2} --validateMappings -o {output}/{"_".join(read.split("_")[:7])} -p 12'
            os.system(cmd)
        elif read.endswith('.1.gz'):
            read_1 = f'{reads_folder}/{read}'
            read_2 = f'{reads_folder}/{read.replace(".1.gz", ".2.gz")}'
            cmd = f'salmon quant -i {salmon_index} -l A -1 {read_1} -2 {read_2} --validateMappings -o {output}/{"_".join(read.split("_")[:7])} -p 12'
            os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: salmon_quantifier.py <reads_folder> <salmon_index> <outdir>')
    else:
        quantify_transcripts(reads_folder=sys.argv[1], salmon_index=sys.argv[2], output=sys.argv[3])

