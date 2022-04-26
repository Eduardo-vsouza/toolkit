import os


def download():
    # start = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR13687743/SRR13687743'
    start = 'https://sra-download.ncbi.nlm.nih.gov/traces/sra58/SRR/013366/SRR13687744'
    splat = start.split("/")
    name = splat[-1].split("SRR")
    for i in range(26):
        number = int(name[-1])
        number += i
        new = f'{name[0]}{number}'
        full = f'{"/".join(splat[:-1])}/{new}'
        print(full)
        # os.system(f'wget {full}')

download()