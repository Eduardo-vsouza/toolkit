import sys

from Bio import SeqIO


class Filtering(object):
    def __init__(self, fasta, gtf):
        self.fasta = fasta
        self.gtf = gtf

        self.filtered = []

    def get_filtered_entries(self):

        records = SeqIO.parse(self.fasta, 'fasta')
        for record in records:
            self.filtered.append(str(record.description))

    def filter_gtf(self, output):
        new_lines = []
        i = 0
        with open(self.gtf, 'r') as handler, open(output, 'w') as out:
            lines = handler.readlines()
            for line in lines:
                i += 1
                print(f'{i/len(lines)*100}'[:5], end='\r')
                if not line.startswith("#"):
                    cols = line.split("\t")
                    attrs = cols[8].split(";")
                    for a in attrs:
                        if 'orf_id' in a:

                            orf_id = a.split(" ")[1].replace("\"", "")
                            # print(orf_id)
                            if orf_id in self.filtered:
                                new_lines.append(line)
                            break
                else:
                    new_lines.append(line)
            out.writelines(new_lines)

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: filter_gtf_by_fasta.py <fasta> <gtf> <output>")
    else:
        data = Filtering(fasta=sys.argv[1], gtf=sys.argv[2])
        data.get_filtered_entries()
        data.filter_gtf(output=sys.argv[3])
