import pandas as pd
from Bio.Blast import NCBIXML


class BlastParser:
    def __init__(self, xml):
        self.xml = xml

    def parse(self, evalue=0.01, score=50):
        for record in NCBIXML.parse(open(self.xml)):
            if record.alignments:  # skip queries with no matches
                for align in record.alignments:

                    for hsp in align.hsps:
                        if hsp.expect <= evalue and hsp.score >= score:
                            print(vars(align))
                            # if (sbjct_end - sbjct_start) >
                            # print([var for var in vars(hsp)])
                            break
            break


if __name__ == '__main__':
    data = BlastParser(xml='/media/eduardo/gold/Eduardo/TT_cells/ribo-seq/ribocode/tt_smorfs_blasted.xml')
    data.parse()