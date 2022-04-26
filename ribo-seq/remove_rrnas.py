import sys
import os

class RibosomalRemover(object):
    def __init__(self, reads_folder):
        self.readsFolder = reads_folder
        self.readsFiles = os.listdir(reads_folder)

    def remove_ribosomal(self):
        cmd = f'bowtie2 -p 8 '