import sys
import os


class BatchPreFetcher:
    def __init__(self, accession_list):
        self.accessionListFile = accession_list

        self.accessions = []

    def get_accessions(self):
        with open(self.accessionListFile, 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                line = line.rstrip()
                self.accessions.append(f'{line}.sra')

    def download(self):
        for acc in self.accessions:
            cmd = f'prefetch {acc}'
            os.system(cmd)

    def dump(self):
        folders = os.listdir(".")
        for folder in folders:
            if os.path.isdir(folder):
                sra = f'{folder}/{folder}.sra'
                cmd = f'fasterq-dump --outdir fastq --skip-technical ' \
                      f'--split-3 -e 12 {sra}'
                os.system(cmd)


if __name__ == '__main__':
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("usage: .py <accession_list_file>\n"
              "Downloads SRA files from a .txt file containing one accession per line, such as:\n"
              "SRR12345\n"
              "SRR12346\n"
              "and uses fastq dump to automatically generate the compressed .fastq.gz files from those. Be aware that"
              " the files are going to be downloaded to the directory where the script was executed.")
    else:
        data = BatchPreFetcher(accession_list=sys.argv[1])
        data.get_accessions()
        data.download()
        data.dump()
