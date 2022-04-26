#!/usr/bin/env python3

import os
import sys


class CutAdapt(object):
    def __init__(self, folder, outdir):
        self.folder = folder
        self.outdir = outdir
        self.check_dir()
        self.adaptor = 'AGATCGGAAGAG'
        self.reverseAdaptor = 'CTCTTCCGATCT'

    def trim(self):
        files = os.listdir(self.folder)
        for file in files:
            if 'fastq' in file and '_R1_' in file:

                pair = file.replace("_R1_", "_R2_")
                cmd = f'cutadapt -a {self.adaptor} -A {self.adaptor} -m 15 -j 12 -o {self.outdir}/{file} ' \
                      f'-p {self.outdir}/{pair} {self.folder}/{file} {self.folder}/{pair}'
                os.system(cmd)

    def check_dir(self):
        if not os.path.exists(self.outdir):
            os.system(f'mkdir {self.outdir}')

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: cutadaptor.py <folder> <outdir>')
    else:
        data = CutAdapt(folder=sys.argv[1], outdir=sys.argv[2])
        data.trim()
