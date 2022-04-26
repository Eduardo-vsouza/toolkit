#!/usr/bin/env python3

import os
import sys

import pandas as pd


class ChosenSequencesReplacer(object):
    def __init__(self, file, output):
        self.df = pd.read_csv(file, sep='\t')
        self.output = output

    def replace(self):
        self.df = self.df.drop(columns=["chosen_sequences"])
        self.df = self.df.rename(columns={"new_chosen_sequences": "chosen_sequences"})

    def save(self):
        self.df.to_csv(self.output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('Replaces the old "chosen_sequences" for the "new_chosen_sequences", updating its name to '
              '"chosen_sequences". \n'
              'usage: replace_old_chosen_sequences.py <file> <output>')
    else:
        data = ChosenSequencesReplacer(file=sys.argv[1], output=sys.argv[2])
        data.replace()
        data.save()
