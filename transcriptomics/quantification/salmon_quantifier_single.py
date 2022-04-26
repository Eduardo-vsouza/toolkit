import sys
import os


def quantify_transcripts(reads_folder, salmon_index, output):
    """ Aligns short reads illumina data to a transcriptome fasta file (not to the genome) and quantifies it using
    salmon. """
    reads = os.listdir(reads_folder)
    for read in reads:
        # if '_R1_' in read and read.endswith('.fq.gz'):
        if 'fastq' in read:
            # read_1 = f'{reads_folder}/{read}'
            # read_2 = f'{reads_folder}/{read.replace("_R1_", "_R2_").replace("val_1", "val_2")}'
            cmd = f'salmon quant -i {salmon_index} -l A -r {reads_folder}/{read} --validateMappings -o {output}/{read}_quantified -p 12'
            os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: salmon_quantifier.py <reads_folder> <salmon_index> <outdir>')
    else:
        quantify_transcripts(reads_folder=sys.argv[1], salmon_index=sys.argv[2], output=sys.argv[3])

