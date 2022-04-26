import sys

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn2


class PipelineComparison(object):
    def __init__(self, riborf_results, ribocode_results):
        self.riborf = riborf_results
        self.ribocode = ribocode_results

        self.ribocodeViral, self.ribocodeHuman = self.get_smorfs(self.ribocode)
        self.riborfViral, self.riborfHuman = self.get_smorfs(self.riborf)

    def get_smorfs(self, fasta):
        viral = []
        human = []
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            name = str(record.description)
            if 'PCMN' in name:
                viral.append(str(record.seq))
            else:
                human.append(str(record.seq))
        return viral, human

    def draw_venns(self):
        ribocode_human = set(self.ribocodeHuman)
        riborf_human = set(self.riborfHuman)
        # set_421 = set(self.humanORFs['421'])
        # set_ui = set(self.humanORFs['UI'])
        venn2([ribocode_human, riborf_human], ('Ribocode', 'Riborf'))
        plt.show()

    def draw_viral_venns(self):
        ribocode_v = set(self.ribocodeViral)
        riborf_v = set(self.riborfViral)
        venn2([ribocode_v, riborf_v], ('Ribocode', 'RibORF'))
        plt.show()

if __name__ == '__main__':
    data = PipelineComparison(riborf_results=f'/media/eduardo/DATA/Eduardo/adenovirus/ribo-seq/riborf/RiboCodeGTF_RiboSeqPipeline_hg19+PCMN102_corr/210706_SE51_Sag_Cindy_hg19.pep',
                              ribocode_results=f'/media/eduardo/DATA/Eduardo/adenovirus/ribo-seq/ribocode/102_reference/summarized_results/humna_and_viral_novel_smorfs.fasta')
    # data.draw_venns()
    data.draw_viral_venns()