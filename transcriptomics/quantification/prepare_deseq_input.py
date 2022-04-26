import os
import sys

import pandas as pd


class DESeqInput(object):
    def __init__(self, folder, output):
        self.folder = folder
        self.samples = {}
        self.output = output
        self.geneList = []

        self.sampleCounts = {}

    def prepare_short_reads(self):
        lines = [',condition,type\n']
        files = os.listdir(self.folder)
        # for subdir in files:
        #     samples = os.listdir(f'{self.folder}/{subdir}')
        for file in files:
            splat = file.split("_")
            print(splat)
            if not 'control' in file:
                sample = '_'.join(splat[:4])
            else:
                sample = '_'.join(splat[:3])
            group = splat[0]
            self.samples[file] = sample
            lib_type = 'single-end'
            # lines.append(f'{sample},{group},{lib_type}\n')
            lines.append(f'{sample},{"_".join(splat[:3])},{lib_type}\n')

            quant = f'{self.folder}/{file}/quant.sf'
            sample_counts = self.__extract_counts(quant)
            self.sampleCounts[sample] = sample_counts

        with open(f'{self.output}_sample_information.txt', 'w') as out:
            out.writelines(lines)

    def create_count_table(self):
        data = {'gene': self.geneList}
        for sample in self.sampleCounts:
            data[sample] = self.sampleCounts[sample]
        df = pd.DataFrame(data)
        df.to_csv(f'{self.output}_counts.txt', sep='\t', index=False)

    def __extract_counts(self, quant_file):
        df = pd.read_csv(quant_file, sep='\t')
        names = df["Name"].tolist()
        if len(self.geneList) == 0:
            self.geneList = names
        counts = df["NumReads"].tolist()
        int_counts = []
        for count in counts:
            int_counts.append(int(count))
        return int_counts


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: <folder> <output_file>')
    else:
        data = DESeqInput(folder=sys.argv[1], output=sys.argv[2])
        data.prepare_short_reads()
        data.create_count_table()