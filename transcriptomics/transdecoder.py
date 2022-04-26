import sys
import os


class TransDecoder(object):
    def __init__(self, assembly_gtf, genome, outdir):
        self.transFolder = '/home/eduardo/programs/TransDecoder/TransDecoder-TransDecoder-v5.5.0/'
        self.assemblyGTF = assembly_gtf
        self.genome = genome
        self.outdir = outdir

    def generate_transcriptome_fasta(self):
        cmd = f'{self.transFolder}/util/gtf_genome_to_cdna_fasta.pl {self.assemblyGTF} {self.genome} ' \
              f'> {self.outdir}/transcripts.fasta'
        os.system(cmd)

    def gtf_to_gff3(self):
        cmd2 = f'{self.transFolder}/util/gtf_to_alignment_gff3.pl {self.assemblyGTF} > {self.outdir}/transcripts.gff3'
        os.system(cmd2)

    def identify_orfs(self):
        cmd3 = f'{self.transFolder}/TransDecoder.LongOrfs -t {self.outdir}/transcripts.fasta'
        os.system(cmd3)

    def last_step(self):
        cmd = f'{self.transFolder}/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder_dir/longest_orfs.gff3 ' \
              f'{self.outdir}/transcripts.gff3 {self.outdir}/transcripts.fasta > {self.outdir}/transcripts.fasta.transdecoder.genome.gff3'
        os.system(cmd)

    def to_gtf(self):
        cmd = f'gffread -T -o {self.outdir}/transcripts_to_update.gtf {self.outdir}/transcripts.fasta.transdecoder.genome.gff3'
        os.system(cmd)

    def update(self):
        cmd = f'GTFupdate {self.outdir}/transcripts_to_update.gtf > {self.outdir}/transcripts_updated_for_ribocode.gtf'
        os.system(cmd)

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print('usage: transdecoder.py <assembly_gtf> <genome> <outdir> \n\n'
              'Be aware that some files are gonna be generated in the working directory.\n')
    else:
        data = TransDecoder(assembly_gtf=sys.argv[1], genome=sys.argv[2], outdir=sys.argv[3])
        data.generate_transcriptome_fasta()
        data.gtf_to_gff3()
        data.identify_orfs()
        data.last_step()
        data.to_gtf()
        data.update()
