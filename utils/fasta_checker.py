from Bio import SeqIO


def check(fasta, output):
    records = SeqIO.parse(fasta, 'fasta')
    new = []
    for record in records:
        entry = str(record.description)
        if '	' in entry:
            splat = entry.split("	")
            new.append(f'>{splat[0]}\n{splat[1]}')
        else:
            new.append(f'>{entry}\n{str(record.seq)}\n')
    with open(output, 'w') as handler:
        handler.writelines(new)


def fix(fasta, output):
    outfile = []
    with open(fasta, 'r') as handler:
        lines = handler.readlines()
        for line in lines:
            if '>' in line:
                if not line.startswith(">"):
                    line = line.rstrip()
                    splat = line.split(">")
                    print(splat)
                    new_entry = f'>{"".join(splat[1:])}'
                    print(new_entry)
                    seq = splat[0]
                    print(seq)
                    new_line = f'{new_entry}\n{seq}\n'
                    outfile.append(new_line)
                else:
                    outfile.append(line)
            else:
                outfile.append(line)
    with open(output, 'w') as out:
        out.writelines(outfile)




if __name__ == '__main__':
    # check('/media/eduardo/gold/Eduardo/ad5_ATCC_Tyger/annotation_files/hg19ad5.fa',
    #       output='/media/eduardo/gold/Eduardo/ad5_ATCC_Tyger/annotation_files/hg19ad5_fixed.fa')
    fix('/media/eduardo/gold/Eduardo/ad5_ATCC_Tyger/annotation_files/hg19ad5_fixed.fa', '/media/eduardo/gold/Eduardo/ad5_ATCC_Tyger/annotation_files/hg19ad5_fixed_2.fa')