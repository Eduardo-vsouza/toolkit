import os
import sys

import pandas as pd
from Bio import SeqIO
from venn import venn
import matplotlib.pyplot as plt


class MSFraggerPostProcessing:
    def __init__(self, df, tags, annotated_tag, outdir):
        self.df = pd.read_csv(df, sep='\t')
        self.df = self.df[self.df["Protein"].str.contains("HUMAN") == False]
        self.df = self.df[self.df["Mapped Proteins"].str.contains("HUMAN") == False]
        self.df = self.df[self.df["Protein"].str.contains("annotated") == False]
        self.df = self.df[self.df["Mapped Proteins"].str.contains("annotated") == False]

        self.tags = tags
        self.annotatedTag = annotated_tag

        self.outdir = outdir
        self.__check_outdir()

        # self.columns = ['Total Spectral Count', 'Unique Spectral Count',
        #                 'Percent Coverage']
        self.columns = ["Spectral Count"]
        self.tagsDataFrames = {tag: {} for tag in self.tags}

    def __check_outdir(self):
        if not os.path.exists(self.outdir):
            os.system(f'mkdir {self.outdir}')

    def separate_proteins(self):
        proteins = self.df["Protein"].tolist()
        mapped = self.df["Indistinguishable Proteins"].tolist()
        # total_specs = self.df["Total Spectral Count"].tolist()
        # uniq_specs = self.df["Unique Spectral Count"].tolist()
        # coverage = self.df["Percent Coverage"].tolist()

        cols = {col: self.df[col].tolist() for col in self.columns}

        for i, protein in enumerate(proteins):
            if type(mapped[i]) == float:
                all_proteins = []
            else:
                all_proteins = mapped[i].split(", ")
            all_proteins.append(protein)
            # for each tag in the file, separate protein entries that contain this tag into another array,
            # along with other pertinent data
            for tag in self.tags:
                for prot in all_proteins:
                    add = True
                    if tag in prot:
                        # if tag == self.annotatedTag:
                        #     if tag not in protein:

                        if add:
                            if 'Protein' not in self.tagsDataFrames[tag]:
                                self.tagsDataFrames[tag]['Protein'] = []
                                self.tagsDataFrames[tag]['Other mapped proteins'] = []
                            self.tagsDataFrames[tag]['Protein'].append(prot)
                            self.tagsDataFrames[tag]['Other mapped proteins'].append(','.join([p for p in all_proteins if p != prot]))

                            for col in cols:
                                if col not in self.tagsDataFrames[tag]:
                                    self.tagsDataFrames[tag][col] = []
                                self.tagsDataFrames[tag][col].append(cols[col][i])

    def separate_peptides(self):
        peptides = self.df["Peptide"].tolist()
        proteins = self.df["Protein"].tolist()
        mapped = self.df["Mapped Proteins"].tolist()
        # total_specs = self.df["Total Spectral Count"].tolist()
        # uniq_specs = self.df["Unique Spectral Count"].tolist()
        # coverage = self.df["Percent Coverage"].tolist()

        cols = {col: self.df[col].tolist() for col in self.columns}

        for i, protein in enumerate(proteins):
            pep = peptides[i]
            if type(mapped[i]) == float:
                all_proteins = []
            else:
                all_proteins = mapped[i].split(", ")
            all_proteins.append(protein)
            # for each tag in the file, separate protein entries that contain this tag into another array,
            # along with other pertinent data
            for tag in self.tags:
                for prot in all_proteins:
                    add = True
                    if tag in prot:
                        # if tag == self.annotatedTag:
                        #     if tag not in protein:
                        if add:
                            if 'Protein' not in self.tagsDataFrames[tag]:
                                self.tagsDataFrames[tag]['Protein'] = []
                                self.tagsDataFrames[tag]['Other mapped proteins'] = []
                                self.tagsDataFrames[tag]['Peptide'] = []
                            self.tagsDataFrames[tag]['Protein'].append(prot)
                            self.tagsDataFrames[tag]['Peptide'].append(pep)
                            self.tagsDataFrames[tag]['Other mapped proteins'].append(','.join([p for p in all_proteins if p != prot]))
                            for col in cols:
                                if col not in self.tagsDataFrames[tag]:
                                    self.tagsDataFrames[tag][col] = []
                                self.tagsDataFrames[tag][col].append(cols[col][i])

    def save(self):
        for tag in self.tagsDataFrames:
            print(tag)
            data = {}
            for col in self.tagsDataFrames[tag]:
                data[col] = self.tagsDataFrames[tag][col]
            df = pd.DataFrame(data=data)
            if tag != self.annotatedTag:
                df = df[df["Other mapped proteins"].str.contains(self.annotatedTag) == False]
            df.to_csv(f'{self.outdir}/{tag}.xls', sep='\t', index=False)

    def generate_fasta_peptides(self):
        outputs = os.listdir(f'{self.outdir}')
        for file in outputs:
            if file.endswith(".xls"):
                fasta = []

                df = pd.read_csv(f'{self.outdir}/{file}', sep='\t')
                pep_number = {}
                proteins = df["Protein"].tolist()
                peptides = df["Peptide"].tolist()
                i = 0
                for prot, pep in zip(proteins, peptides):
                    name = f'{prot}_pep'
                    if name not in pep_number:
                        pep_number[name] = 0
                    pep_number[name] += 1

                    new_entry = f'{name}_{pep_number[name]}'
                    record = f'>{new_entry}\n{pep}\n'
                    fasta.append(record)
                with open(f'{self.outdir}/{file[:-4]}_peptides.fasta', 'w') as handler:
                    handler.writelines(fasta)


    def generate_fasta(self, database):
        files = os.listdir(self.outdir)
        for file in files:
            if file.endswith('.xls'):
                fasta = []
                checker = []
                unique = []
                df = pd.read_csv(f'{self.outdir}/{file}', sep='\t')
                proteins = df["Protein"].tolist()
                records = SeqIO.parse(database, 'fasta')
                for record in records:
                    entry = str(record.description)

                    if entry in proteins:
                        if str(record.seq) not in checker:
                            checker.append(str(record.seq))
                            unique.append(f'>{entry}\n{str(record.seq)}\n')
                        fasta.append(f'>{entry}\n{str(record.seq)}\n')
                with open(f'{self.outdir}/{file[:-4]}.fasta', 'w') as handler, open(f'{self.outdir}/{file[:-4]}_unique.fasta', 'w') as uniq_out:
                    handler.writelines(fasta)
                    uniq_out.writelines(unique)

    def create_venns(self):
        files = os.listdir(self.outdir)
        prots = {}
        for file in files:
            if file.endswith("unique.fasta"):
                name = file.split(".")[0].replace("_unique", "")
                if name not in prots:
                    prots[name] = set()
                records = SeqIO.parse(f'{self.outdir}/{file}', 'fasta')
                for record in records:
                    prots[name].add(str(record.seq))
        venn(prots)
        # plt.show()
        plt.savefig(f'{self.outdir}/venn_diagram.png')



if __name__ == '__main__':
    # folder = '/media/eduardo/gold/Eduardo/adenovirus/new_replicates/102_421_tyger_nanopores_matthews_human_uniprot_all_3ft_only_virus/1_1/'
    # folder = '/media/eduardo/gold/Eduardo/adenovirus/new_replicates/102_421_tyger_nanopores_matthews_human_uniprot_all_3ft_only_virus/2_2/'
    folder = f'/media/eduardo/gold/Eduardo/adenovirus/new_replicates/25_04_22_search_102_421_tyger_nanopores_matthews_human_uniprot_all_3ft_only_virus'

    data_102 = MSFraggerPostProcessing(df=f'{folder}/peptide.tsv',
                                       tags=('matthews', '421_nanopore', '102_nanopore', 'tyger'),
                                       annotated_tag='annotated', outdir=f'{folder}/summarized_results')
    # data_102.separate_proteins()

    # data_102.save()
    # data_102.generate_fasta(database='/media/eduardo/gold/Eduardo/adenovirus/proteogenomics/3ft_22_04_22/2022-04-22-decoys-contam-102_421_tyger_nanopore_assemblies_matthews_human_uniprot_all_3ft.fasta.fas')
    # data_102.create_venns()
    data_102.separate_peptides()
    data_102.save()
    data_102.generate_fasta_peptides()

