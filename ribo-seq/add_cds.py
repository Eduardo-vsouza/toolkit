import sys
import os


class CDSPicking(object):
    def __init__(self, gtf_with_cds, gtf_to_add_cds):
        self.gtfWithCDS = gtf_with_cds
        self.targetGTF = gtf_to_add_cds
        self.cdsByGene = {}  # dictionary with key, value: gene_id, [line with a CDs of this gene, second line, 3rd...]

    def parse(self):
        with open(self.gtfWithCDS, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                if not line.startswith("#"):
                    feature, attrs = self.__get_info(line)
                    if feature == 'CDS':
                        for a in attrs:
                            if 'gene_id' in a or 'transcript_id' in a or 'ref_gene_id' in a:
                                gene_id = a.split(" ")[1].replace("\"", "")
                                if gene_id not in self.cdsByGene:
                                    self.cdsByGene[gene_id] = []
                                if line not in self.cdsByGene[gene_id]:
                                    self.cdsByGene[gene_id].append(line)

    def add_cds(self, output):
        new_lines = []
        with open(self.targetGTF, 'r') as handler, open(output, 'w') as outfile:
            lines = handler.readlines()
            for line in lines:
                if not line.startswith("#"):
                    feature, attrs = self.__get_info(line)
                    for a in attrs:
                        if 'gene_id' in a or 'transcript_id' in a or 'ref_gene_id' in a:
                            gene_id = a.split(" ")[1].replace("\"", "")
                            if gene_id in self.cdsByGene:
                                cds_lines = self.cdsByGene[gene_id]
                                for cline in cds_lines:
                                    if cline not in new_lines:
                                        new_lines.append(cline)
                new_lines.append(line)
            outfile.writelines(new_lines)

    @staticmethod
    def __get_info(line):
        cols = line.split("\t")
        feature = cols[2]
        attrs = cols[8].split(";")
        return feature, attrs


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: add_cds.py <gtf_with_cds> <gtf_to_add_cds_to> <output>")
    else:
        data = CDSPicking(gtf_with_cds=sys.argv[1], gtf_to_add_cds=sys.argv[2])
        data.parse()
        data.add_cds(output=sys.argv[3])
