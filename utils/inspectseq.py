from Bio import SeqIO


def check(fasta):
    records = SeqIO.parse(fasta, 'fasta')
    aa = 'ARNDCEQGHOILKMFPUSTWYV'
    aas = [i for i in aa]
    print(aas)
    for record in records:
        # print(str(record.seq))
        print(record.id)
        # print(len(str(record.seq)))
        # print('\n')
        # for i in str(record.seq):
        #     if i not in aas:
        #         print(record.id)
                # print(str(record.seq))
                # print(i)


if __name__ == '__main__':
    check('/media/eduardo/DATA/Eduardo/adenovirus/updated_annotation_files/virusproteome')