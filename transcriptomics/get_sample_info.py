#!/usr/bin/python3

import sys
import os

import pandas as pd


class SampleDataCollector:
    def __init__(self, df):
        self.df = pd.read_csv(df, sep=',')

        self.data = {'run': [], 'condition': []}

    def reformat(self):
        runs = self.df["Run"].tolist()
        group = self.df["drug_treatment"].tolist()
        times = self.df["treatment_time"].tolist()

        for i in range(len(runs)):
            self.data['run'].append(runs[i])
            time = times[i].split(" ")[0]
            condition = f'{group[i].replace(" ", "_").replace("/", "-")}_{time}h'
            self.data['condition'].append(condition)

    def save(self, output):
        df = pd.DataFrame(data=self.data)
        df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print(f"usage: {os.path.basename(__file__)} <df> <output>\n\n"
              f"{os.path.basename(__file__)} takes a data frame <df> from NCBI SRA Run Selector, named "
              f"SraRunTable.txt by default, and generates another file <output> containing the two columns: 'run' and "
              f"'condition'. These can be used to rename the read or quantification files later with "
              f"sample_renamer_quantification.py.")
    else:
        data = SampleDataCollector(df=sys.argv[1])
        data.reformat()
        data.save(output=sys.argv[2])