import os
import sys


class ViralTranscriptome(object):
    def __init__(self, gtf):
        self.gtf = gtf

        self.transcripts = []

    def get_viral(self):
        with open(self.gtf, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if not line.startswith("#"):
                    cols = line.split("\t")
                    chrom = cols[0]
                    if 'PCMN' in chrom:
                        attrs = cols[8].split(";")
                        for a in attrs:
                            if 'transcript_id' in a:
                                name = a.split(" ")[2].replace("\"", "").replace(";", "")
                                if f'{name}\n' not in self.transcripts:
                                    self.transcripts.append(f'{name}\n')

    def save(self, output):
        with open(output, 'w') as handler:
            handler.writelines(self.transcripts)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: viral_transcripts.py <gtf> <output_file')
    else:
        data = ViralTranscriptome(gtf=sys.argv[1])
        data.get_viral()
        data.save(output=sys.argv[2])