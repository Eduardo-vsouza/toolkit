import os
import sys


class RibosomalRemover(object):
    def __init__(self, reads_folder, rrnas, outdir, threads=12, experiment='paired'):
        self.readsFolder = reads_folder
        self.rRNAs = rrnas
        self.threads = threads
        self.outdir = outdir
        self.experiment = experiment

    def index(self):
        cmd = f'bowtie2-build {self.rRNAs} rRNA'
        os.system(cmd)

    def align(self):
        self.__check_dir(self.outdir)
        contamined = f'{self.outdir}/ribosome_contaminated_reads'
        cleaned = f'{self.outdir}/ribosome_cleaned_reads'
        self.__check_dir(contamined)
        self.__check_dir(cleaned)
        reads = os.listdir(self.readsFolder)
        for file in reads:
            if '.fastq' in file:
                file_path = f'{self.readsFolder}/{file}'
                if self.experiment == 'single':
                    cmd = f'bowtie2 -p {self.threads} --norc --un {cleaned}/{file} -U ' \
                          f'{file_path} -x rRNA -S {contamined}/{file}.sam'
                    os.system(cmd)

                elif self.experiment == 'paired':
                    if '_1.' in file:
                        pair = f'{self.readsFolder}/{file.replace("_1.", "_2.")}'
                        cmd = f'bowtie2 -p {self.threads} --norc --un-conc {cleaned}/{file} -1 ' \
                              f'{file_path} -2 {pair} -x rRNA -S {contamined}/{file}.sam'
                        os.system(cmd)
    @staticmethod
    def __check_dir(folder):
        if not os.path.exists(folder):
            os.system(f'mkdir {folder}')


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: ribosomal_remover.py <reads_folder> <rrnas_fasta> <outdir> <threads>')
    else:
        data = RibosomalRemover(reads_folder=sys.argv[1], rrnas=sys.argv[2], outdir=sys.argv[3],
                                threads=int(sys.argv[4]))
        data.index()
        data.align()
