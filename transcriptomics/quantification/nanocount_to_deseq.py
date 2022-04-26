import sys
import os

import pandas as pd


class NanocountToDESeq2(object):
    def __init__(self, nanocounts_folder, outdir):
        self.nanoFolder = nanocounts_folder
        self.outdir = outdir
        self.__check_dir()

        self.transcripts = []
        self.added = False
        self.counts = {}

    def __check_dir(self):
        if not os.path.exists(self.outdir):
            os.system(f'mkdir {self.outdir}')

    def sample_information(self):
        files = os.listdir(self.nanoFolder)
        to_write = [',condition,type\n']
        for file in files:
            splat = file.split("_")
            group = splat[0]
            replicate = '_'.join(splat[:2])
            to_write.append(f'{replicate},{group},single-end\n')
            self.__extract_counts(f'{self.nanoFolder}/{file}', replicate)
        with open(f'{self.outdir}/sample_information.txt', 'w') as out:
            out.writelines(to_write)

    def __extract_counts(self, file, rep):
        df = pd.read_csv(file, sep='\t')
        # if rep not in self.counts:
            # self.counts[rep] = []
        counts = df["est_count"].tolist()
        transcripts = df["transcript_name"].tolist()
        if not self.added:
            for transcript in transcripts:
                if transcript not in self.counts:
                    self.counts[transcript] = {rep: []}
        self.added = True

        # for transcript in self.counts:
        #     if rep not in self.counts[transcript]:
        #         self.counts[transcript][rep] =
        for transcript, count in zip(transcripts, counts):
            if transcript in self.counts:
                self.counts[transcript][rep] = int(count)

        for transcript in self.counts:
            if rep not in self.counts[transcript]:
                self.counts[transcript][rep] = 0


    def counts_table(self):
        data = {'gene': []}
        for transcript in self.counts:
            data['gene'].append(transcript)
            for rep in self.counts[transcript]:
                if rep not in data:
                    data[rep] = []
                data[rep].append(self.counts[transcript][rep])


        df = pd.DataFrame(data=data)
        df.to_csv(f'{self.outdir}/counts_table.txt', sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: nanocount_to_deseq.py <nanocount_folder> <output_directory>')
    else:
        data = NanocountToDESeq2(nanocounts_folder=sys.argv[1], outdir=sys.argv[2])
        data.sample_information()
        data.counts_table()
