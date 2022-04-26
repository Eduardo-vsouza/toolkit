import sys

import pandas as pd


class GTFGenerator:
    def __init__(self, df, output, tiers=('T1', 'T2', 'T3')):
        self.df = pd.read_csv(df, sep='\t')
        self.output = output
        self.tiers = tiers

    def create_gtf(self):
        entries = self.df["Final Entries"].tolist()
        coords = self.df["Genome Coordinates"].tolist()
        tiers = self.df["Tier"].tolist()
        to_write = []
        checker = []
        # <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        for i, entry in enumerate(entries):
            if tiers[i] in self.tiers:
                if entry not in checker:
                    checker.append(entry)
                    chromosome = 'NC_000962.3'
                    source = 'uproteins'
                    score = '.'
                    frame = 0
                    feature = 'CDS'
                    start, end, strand = self.__define_coordinates(coords[i], entry)
                    gene_id = entry
                    attributes = f'gene_id \"{gene_id}\"; locus_tag={"_".join(entry.split("_")[:3])}; ' \
                                 f'protein_id={"_".join(entry.split("_")[:3])}'
                    line = f'{chromosome}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}\n'
                    to_write.append(line)
        with open(self.output, 'w') as handler:
            handler.writelines(to_write)

    @staticmethod
    def __define_coordinates(coordinates, entry):
        splat_coords = coordinates.split("-")
        if 'forward' in entry:
            start = int(splat_coords[0])
            end = int(splat_coords[1])
            strand = '+'
        elif 'reverse' in entry:
            start = int(splat_coords[0])
            end = int(splat_coords[1])
            strand = '-'
        else:
            start = int(splat_coords[0])
            end = int(splat_coords[1])
            strand = '+'
            if start > end:
                end, start = start, end
                strand = '-'
        return start, end, strand


if __name__ == '__main__':
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("usage: .py <df> <output_gtf>")
    else:
        data = GTFGenerator(df=sys.argv[1], output=sys.argv[2])
        data.create_gtf()
