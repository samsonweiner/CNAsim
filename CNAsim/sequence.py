from pyfaidx import Fasta
import numpy as np
import subprocess
import os
import pickle

def read_fasta(input_fasta, chrom_lens_only = False):
    ref = Fasta(input_fasta, one_based_attributes=False)
    chrom_lens = {}
    for chrom in ref.keys():
        chrom_lens[chrom] = len(ref[chrom])
    if chrom_lens_only:
        return chrom_lens
    else:
        return ref, chrom_lens

def make_alt_ref(ref):
    pass

def build_cell_ref(genome, ref, chrom_names, regions_per_chrom, region_length, allele, prefix):
    if isinstance(genome, str):
        with open(genome, 'rb') as f:
            genome = pickle.load(f)
    ref_name = prefix + '_allele' + str(allele) + '.fa'
    chrom_lens = {}
    f = open(ref_name, 'w+')
    for chrom in chrom_names:
        chrom_lens[chrom] = 0
        num_regions = regions_per_chrom[chrom]
        if len(genome[chrom][allele]) > 0:
            f.write('>' + chrom + '\n')
            for homolog in genome[chrom][allele]:
                while 'X' in homolog:
                    homolog.remove('X')
                for r in homolog:
                    if r == num_regions - 1:
                        f.write(ref[chrom][r*region_length:].seq)
                        chrom_lens[chrom] += len(ref[chrom][r*region_length:].seq)
                    else:
                        f.write(ref[chrom][r*region_length:(r+1)*region_length].seq)
                        chrom_lens[chrom] += region_length
            f.write('\n')
    f.close()
    call = subprocess.run(['samtools', 'faidx', ref_name])
    #cur_ref = Fasta(ref_name)
    return ref_name, chrom_lens
