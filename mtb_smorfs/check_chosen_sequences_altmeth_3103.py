import sys

import pandas as pd


def compare(old, new, output):
    old_df = pd.read_csv(old, sep='\t')
    new_df = pd.read_csv(new, sep='\t')
    old_chosen = old_df["chosen_sequences"].tolist()
    new_df = new_df[new_df["chosen_sequences"].isin(old_chosen)]
    new_df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    compare(old=sys.argv[1], new=sys.argv[2], output=sys.argv[3])