import sys


class RibocodeGTFFormatter(object):
    def __init__(self, gtf):
        self.gtf = gtf

    def parse(self, output):
        new_lines = []
        with open(self.gtf, 'r') as handler, open(output, 'w') as outfile:
            lines = handler.readlines()
            for line in lines:
                if line.startswith('PCMN'):
                    if 'gene_id' not in line:
                        line = line.rstrip()
                        cols = line.split("\t")
                        attrs = cols[8].split(";")
                        new_attrs_list = []
                        for a in attrs:
                            if 'transcript_id' in a:
                                gene = a.split(" ")[1]
                                new_attrs_list.insert(0, f' gene_id \"{gene}\"')
                                for i in attrs:
                                    if 'gene_id' not in i:
                                        new_attrs_list.append(i)
                                # new_attrs_list.append(f' gene_id \"{gene}\";\n')

                        new_attrs = ';'.join(new_attrs_list).replace("transcript_id", " transcript_id")
                        new_attrs += '\n'
                        old_cols = cols[:-1]
                        old_cols.append(new_attrs)
                        new_cols = '\t'.join(old_cols)
                        new_lines.append(new_cols)
                    else:
                        new_lines.append(line)
                else:
                    new_lines.append(line)
            outfile.writelines(new_lines)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: ribocode_gtf.py <gtf> <output>")
    else:
        data = RibocodeGTFFormatter(gtf=sys.argv[1])
        data.parse(output=sys.argv[2])

