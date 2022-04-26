import pandas as pd


def get_ribo(tsv):
    df = pd.read_csv(tsv, sep='\t')
    df = df[df["Feature"] == 'rRNA']
    print(df["Start"])
    print(df["Stop"])
    print(df["Strand"])
    print(df["Product"])

get_ribo('/home/eduardo/Downloads/Mycobacterium_tuberculosis_H37Rv_txt_v4.txt')