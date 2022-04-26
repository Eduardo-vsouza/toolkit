import sys
import os


class TrimmomaticIterator(object):
    def __init__(self, folder, outdir):
        self.folder = folder
        self.outdir = outdir
        self.__check_dir(self.folder)
        self.paired = f'{self.outdir}/paired'
        self.unpaired = f'{self.outdir}/unpaired'
        self.__check_dir(self.paired)
        self.__check_dir(self.unpaired)

    @staticmethod
    def __check_dir(folder):
        if not os.path.exists(folder):
            os.system(f'mkdir {folder}')

    def trim(self):
        files = os.listdir(self.folder)
        for file in files:
            if '_R1_' in file and 'fastq' in file:
                rep1 = file
                rep2 = file.replace("_R1_", "_R2_")
                cmd = f'java -jar /home/eduardo/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 12 ' \
                      f'{self.folder}/{rep1} {self.folder}/{rep2} {self.paired}/{rep1} {self.unpaired}/{rep1} {self.paired}/{rep2} ' \
                      f'{self.unpaired}/{rep2} ILLUMINACLIP:/home/eduardo/programs/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36'
                os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: trimmomatic.py <folder> <outdir>')
    else:
        data = TrimmomaticIterator(folder=sys.argv[1], outdir=sys.argv[2])
        data.trim()
