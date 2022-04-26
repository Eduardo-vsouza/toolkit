import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class LengthTPMCorrelation(object):
    def __init__(self, gtf_list_comma_sep):
        self.gtfFiles = gtf_list_comma_sep.split(",")
        self.geneLengths = []
        self.tpms = []

    def parse_gtf(self):

        for file in self.gtfFiles:

            with open(file, 'r') as handler:
                lines = handler.readlines()
                for line in lines:
                    if not line.startswith("#"):
                        cols = line.split("\t")
                        if cols[2] == 'transcript':
                            attrs = cols[8].split(";")

                            start = int(cols[3])
                            end = int(cols[4])
                            length = (end - start)
                            # print(length)
                            # if length < 1000:
                            for at in attrs:
                                if at.startswith(" TPM"):
                                    tpm = (float(at.split(" ")[2].replace("\"", "")))
                                    self.tpms.append(tpm)
                            self.geneLengths.append(length)

    def frame(self):
        self.tpms = np.log(self.tpms)
        self.geneLengths = np.log(self.geneLengths)
        self.df = pd.DataFrame(data={'tpms': self.tpms, 'transcript_length': self.geneLengths})


    def plot(self):
        # print(self.geneLengths)
        # sns.histplot(x='transcript_length', y='tpms', data=self.df)
        # sns.displot(x='transcript_length', y='tpms', data=self.df)

        sns.histplot(data=self.df, x='transcript_length')
        plt.show()


if __name__ == '__main__':
    # folder = '/media/eduardo/DATA/Eduardo/adenovirus/nanopore/corrected_transcripts_assemblies/stringtie/102/102_1/hybrid_assembly/'
    folder = '/media/eduardo/DATA/Eduardo/adenovirus/nanopore/corrected_transcripts_assemblies/stringtie/421/421_1/hybrid_assembly'
    # data = LengthTPMCorrelation(gtf_list_comma_sep=f'{folder}/hybrid_assembly_210624_102_1_mRNA_S3.adaptQualTrim._STARalign_hg19Gencode+PCMN102_sort.gtf,{folder}/hybrid_assembly_210624_102_2_mRNA_S4.adaptQualTrim._STARalign_hg19Gencode+PCMN102_sort.gtf')
    data = LengthTPMCorrelation(gtf_list_comma_sep=f'{folder}/hybrid_assembly_210624_421_1_mRNA_S5.adaptQualTrim._STARalign_hg19Gencode+PCMN421_sort.gtf,{folder}/hybrid_assembly_210624_421_2_mRNA_S6.adaptQualTrim._STARalign_hg19Gencode+PCMN421_sort.gtf')
    data.parse_gtf()
    data.frame()
    data.plot()
