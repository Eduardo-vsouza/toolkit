import os
import sys


class HumanRemover(object):
    def __init__(self, reads_folder, index, outdir, threads=12, experiment='paired'):
        self.readsFolder = reads_folder
        self.threads = threads
        self.outdir = outdir
        self.experiment = experiment
        self.index = index

    # def index(self):
    #     cmd = f'hisat2-build {self.rRNAs} rRNA'
    #     os.system(cmd)

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
                        print(f'Removing human reads from {file}')
                        pair = f'{self.readsFolder}/{file.replace(".1.", ".2.")}'
                        cmd = f'hisat2 -p {self.threads} --norc --un-conc {cleaned}/{file} -1 ' \
                              f'{file_path} -2 {pair} -x {self.index} -S {contamined}/{file}.sam'
                        os.system(cmd)
                        print(f'Removed human reads form {file}')

    @staticmethod
    def __check_dir(folder):
        if not os.path.exists(folder):
            os.system(f'mkdir {folder}')


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: ribosomal_remover.py <reads_folder> <index> <outdir> <threads>')
    else:
        data = HumanRemover(reads_folder=sys.argv[1], index=sys.argv[2], outdir=sys.argv[3],
                                threads=int(sys.argv[4]))
        data.align()
