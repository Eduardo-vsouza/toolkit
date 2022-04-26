import sys

import pandas as pd


class FilterThemOut(object):
    def __init__(self, df, output):
        self.df = pd.read_csv(df, sep='\t')
        self.output = output

    def filter(self):
        new_entries = []
        altmeth_seqs = {}
        entries = self.df["Final Entries"].tolist()
        altmeth = self.df["alternative_methionine_sequence"].tolist()
        seqs = self.df["chosen_sequences"].tolist()
        for i in range(len(entries)):
            if seqs[i] == altmeth[i]:
                altmeth_seqs[seqs[i]] = entries[i]
            # else:
            #     altmeth[seqs[i]] =

        for i in range(len(entries)):
            if altmeth[i] != 'False':
                new_entries.append(altmeth_seqs[altmeth[i]])
            else:
                new_entries.append(entries[i])
        df = self.df.drop(columns=['Final Entries'])
        df.insert(3, "Final Entries", new_entries)
        df.to_csv(self.output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("usage: refilter_altmeth.py <df> <output>")
    else:
        data = FilterThemOut(df=sys.argv[1], output=sys.argv[2])
        data.filter()