from pyfaidx import Fasta
import numpy as np
import subprocess

def read_fasta(input_fasta):
    ref = Fasta(input_fasta)
    chrom_lens = {}
    for chrom in ref.keys():
        chrom_lens[chrom] = len(ref[chrom])
    return ref, chrom_lens

def build_cell_ref(genome, ref, min_cn_len, allele, prefix):
    ref_name = prefix + '_allele' + str(allele) + '.fa'
    chrom_lens = {}
    f = open(ref_name, 'w+')
    for chrom in ref.keys():
        chrom_lens[chrom] = 0
        num_regions = round(len(ref[chrom])/min_cn_len)
        f.write('>' + chrom + '\n')

        for r in genome[chrom][allele]:
            if r == num_regions - 1:
                f.write(ref[chrom][r*min_cn_len:].seq)
                chrom_lens[chrom] += len(ref[chrom][r*min_cn_len:].seq)
            else:
                f.write(ref[chrom][r*min_cn_len:(r+1)*min_cn_len].seq)
                chrom_lens[chrom] += min_cn_len
        f.write('\n')
    f.close()
    call = subprocess.run(['samtools', 'faidx', ref_name])
    #ref, chrom_lens = read_fasta(ref_name)
    return ref_name, chrom_lens
