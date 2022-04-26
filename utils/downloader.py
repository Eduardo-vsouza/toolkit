import os
import sys


def download(file):
    with open(file, 'r') as handler:
        lines = handler.readlines()
        for line in lines:
            line = line.rstrip()
            os.system(f'wget {line}')


if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: downloader.py <file>')
    else:
        download(file=sys.argv[1])