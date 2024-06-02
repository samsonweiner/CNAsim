import os
import subprocess
import itertools
from collections import Counter
import pickle
from itertools import repeat
import multiprocessing
from pyfaidx import Fasta
from Bio import SeqIO

import numpy as np
from scipy.stats import beta, poisson
from scipy.optimize import newton_krylov
from scipy.optimize.nonlin import NoConvergence

from .sequence import *

def iter_by_chunk(iterable, chunksize):
    return itertools.zip_longest(*[iter(iterable)] * chunksize)

def init_pool_processes():
    np.random.seed()

# Estimating alpha/beta from point on lorenz curve
def get_alpha_beta(x0, y0):
    def F(P):
        Aa, Bb = P[0], P[1]
        X = beta.cdf(float(Aa)/(Aa+Bb), Aa, Bb) - x0
        Y = beta.cdf(float(Aa)/(Aa+Bb), Aa+1, Bb) - y0
        return [X, Y]
    guess = [10, 5]
    max_n = 1000
    n = 1
    while n < max_n:
        try:
            sol = newton_krylov(F, guess, method = 'lgmres', verbose = 0, rdiff = 0.1, maxiter=50)
            break
        except NoConvergence as e:
            guess = np.random.rand(2) * 10 + 0.1
        except ValueError as e:
            guess = np.random.rand(2) * 10 + 0.1
        n += 1
    if n == max_n:
        # Error has occurred here
        return [1.38, 1.38]
    else:
        return sol

def bezier_coef(points):
    n = len(points) - 1

    C = 4 * np.identity(n)
    np.fill_diagonal(C[1:], 1)
    np.fill_diagonal(C[:, 1:], 1)
    C[0, 0] = 2
    C[n - 1, n - 1] = 7
    C[n - 1, n - 2] = 2

    P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
    P[0] = points[0] + 2 * points[1]
    P[n - 1] = 8 * points[n - 1] + points[n]

    A = np.linalg.solve(C, P)
    B = [0] * n
    for i in range(n - 1):
        B[i] = 2 * points[i + 1] - A[i + 1]
    B[n - 1] = (A[n - 1] + points[n]) / 2

    return A, B

def single_cubic_bezier(a, b, c, d):
    return lambda t: np.power(1 - t, 3) * a + 3 * np.power(1 - t, 2) * t * b + 3 * (1 - t) * np.power(t, 2) * c + np.power(t, 3) * d

def multiple_cubic_bezier(points):
    A, B = bezier_coef(points)
    return [single_cubic_bezier(points[i], A[i], B[i], points[i + 1]) for i in range(len(points) - 1)]
    
# Generates starting coverage points at set intervals 
def gen_start_interval(num_windows, interval, Aa, Bb):
    # Sample point every interval bin
    x = [i for i in range(0, num_windows, interval)]
    if x[-1] != num_windows-1:
        if num_windows - x[-1] >= interval/2:
            x.append(num_windows-1)
        else:
            x[-1] = num_windows-1
    # draw coverage from beta distribution
    y = list(np.random.beta(Aa, Bb, len(x)))
    points = [np.array(p) for p in zip(x, y)]
    return points

# Smooths points across bins with bezier curves
def gen_coverage(num_windows, interval, Aa, Bb):
    points = gen_start_interval(num_windows, interval, Aa, Bb)
    A, B = bezier_coef(points)
    curves = [single_cubic_bezier(points[i], A[i], B[i], points[i + 1]) for i in range(len(points) - 1)]

    new_points = []
    for i in range(len(points) - 1):
        f = curves[i]
        gaps = points[i+1][0] - points[i][0] + 1
        coords = [f(t) for t in np.linspace(0, 1, int(gaps))]
        if i == 0:
            new_points.append((round(coords[0][0]), coords[0][1]))
        new_points += [(round(i[0]), i[1]) for i in coords[1:]]
    
    new_points.sort(key=lambda x: x[0])
    # Multiply by 2 so that values are in range 0-2. Value of 1 indicates mean coverage.
    new_points = [2*max(min(i[1], 1), 0) for i in new_points]
    return new_points

# Draw readcounts for each bin
def draw_readcounts(num_windows, window_size, interval, Aa, Bb, coverage, read_len):
    avg_read = (coverage*window_size) / (2*read_len)
    cov_scales = gen_coverage(num_windows, interval, Aa, Bb)
    #exp_counts = [avg_read*x for x in cov_scales]
    #readcounts = [np.random.poisson(x) for x in exp_counts]
    readcounts = [poisson.rvs(avg_read*x) for x in cov_scales]
    return readcounts

# Generate reads across the genome for a given cell
def gen_reads_cell(cell, chrom_names, uniform_coverage, window_size, interval, Aa, Bb, coverage, read_len, seq_error, out_path):
    #print('Generating reads for:', cell.name)
    prefix = os.path.join(out_path, cell.name)
    cell_ref1 = prefix + '_allele0.fa'
    cell_ref2 = prefix + '_allele1.fa'
    counts = []
    for allele, cell_ref in enumerate([cell_ref1, cell_ref2]):
        init = False
        chrom_lens = read_fasta(cell_ref, chrom_lens_only = True)
        for chrom in chrom_names:
            if chrom in chrom_lens:
                num_windows = round(chrom_lens[chrom] / window_size)
                if not uniform_coverage:
                    readcounts = draw_readcounts(num_windows, window_size, interval, Aa, Bb, coverage, read_len)
                for w in range(num_windows):
                    start = w * window_size + 1
                    if w == num_windows - 1:
                        end = chrom_lens[chrom]
                    else:
                        end = (w+1) * window_size
                    
                    region_fa_path = cell_ref[:-3] + chrom + '-' + str(start) + '.fa'
                    with open(region_fa_path, 'w+') as f:
                        call = subprocess.run(['samtools', 'faidx', cell_ref, chrom + ':' + str(start) + '-' + str(end)], stdout=f)
                    
                    for record in SeqIO.parse(region_fa_path, 'fasta'):
                        total_N = record.seq.count('N')
                    total_bp = end - start + 1
                    N_ratio = total_N / total_bp
                    if uniform_coverage:
                        readcount = round(poisson.rvs((total_bp*coverage) / (2*read_len)) * (1-N_ratio))
                    else:
                        readcount = round(readcounts[w] * (1-N_ratio))

                    counts.append((allele, chrom, start, end, readcount))

                    if readcount > 0:
                        #seed = np.random.randint(1, 10000000)
                        proc = subprocess.run(['dwgsim', '-H', '-o', '1', '-N', str(readcount), '-1', str(read_len), '-2', str(read_len), '-e', str(seq_error), '-E', str(seq_error), region_fa_path, region_fa_path[:-3]], capture_output=True, text=True)
                        if not init:
                            os.rename(region_fa_path[:-3] + '.bwa.read1.fastq.gz', cell_ref[:-3] + '.read1.fastq.gz')
                            os.rename(region_fa_path[:-3] + '.bwa.read2.fastq.gz', cell_ref[:-3] + '.read2.fastq.gz')
                            init = True
                        else:
                            with open(cell_ref[:-3] + '.read1.fastq.gz', 'a') as f1, open(cell_ref[:-3] + '.read2.fastq.gz', 'a') as f2:
                                call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read1.fastq.gz'], stdout=f1)
                                call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read2.fastq.gz'], stdout=f2)
                            os.remove(region_fa_path[:-3] + '.bwa.read1.fastq.gz')
                            os.remove(region_fa_path[:-3] + '.bwa.read2.fastq.gz')
                        os.remove(region_fa_path[:-3] + '.mutations.txt')
                        os.remove(region_fa_path[:-3] + '.mutations.vcf')
                    os.remove(region_fa_path)
                    #os.remove(region_fa_path + '.fai')
        os.remove(cell_ref)
        os.remove(cell_ref + '.fai')
    with open(prefix + '.read1.fastq.gz', 'w+') as f1, open(prefix + '.read2.fastq.gz', 'w+') as f2:
        call = subprocess.run(['cat', prefix + '_allele0.read1.fastq.gz', prefix + '_allele1.read1.fastq.gz'], stdout=f1)
        call = subprocess.run(['cat', prefix + '_allele0.read2.fastq.gz', prefix + '_allele1.read2.fastq.gz'], stdout=f2)
    os.remove(prefix + '_allele0.read1.fastq.gz')
    os.remove(prefix + '_allele1.read1.fastq.gz')
    os.remove(prefix + '_allele0.read2.fastq.gz')
    os.remove(prefix + '_allele1.read2.fastq.gz')
    os.remove(prefix + '.pkl')

    with open(os.path.join(out_path, cell.name + '.readcounts.tsv'), 'w+') as f:
        for c in counts:
            f.write(f'{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\t{c[4]}\n')

# Generate reads for all cells
def gen_reads(ref1, ref2, num_regions, chrom_names, tree, uniform_coverage, x0, y0, region_length, interval, window_size, coverage, read_len, seq_error, out_path, num_processors):
    leaves = list(tree.iter_leaves())

    [Aa, Bb] = get_alpha_beta(x0, y0)
    if num_processors == 1:
        for cell in leaves:
            prefix = os.path.join(out_path, cell.name)
            for allele, cur_ref in enumerate([ref1, ref2]):
                cell_ref, chrom_lens = build_cell_ref(prefix + '.pkl', cur_ref, chrom_names, num_regions, region_length, allele, prefix)
            gen_reads_cell(cell, chrom_names, uniform_coverage, window_size, interval, Aa, Bb, coverage, read_len, seq_error, out_path)
    else:
        pool = multiprocessing.Pool(processes=num_processors)
        for chunk in iter_by_chunk(tree.iter_leaves(), num_processors):
            cells = list(chunk)
            while None in cells:
                cells.remove(None)
            for cell in cells:
                prefix = os.path.join(out_path, cell.name)
                for allele, cur_ref in enumerate([ref1, ref2]):
                    cell_ref, chrom_lens = build_cell_ref(prefix + '.pkl', cur_ref, chrom_names, num_regions, region_length, allele, prefix)
            with multiprocessing.Pool(processes=num_processors) as pool:
                pool.starmap(gen_reads_cell, zip(cells, repeat(chrom_names), repeat(uniform_coverage), repeat(window_size), repeat(interval), repeat(Aa), repeat(Bb), repeat(coverage), repeat(read_len), repeat(seq_error), repeat(out_path)))

def gen_readcounts(tree, chrom_names, bins, num_regions, region_length, uniform_coverage, x0, y0, interval, window_size, coverage, read_len, out_path):
    regions_per_window = round(window_size / region_length)
    [Aa, Bb] = get_alpha_beta(x0, y0)
    leaves = list(tree.iter_leaves())
    final_readcounts = {}

    for cell in leaves:
        prefix = os.path.join(out_path, cell.name)
        with open(prefix + '.pkl', 'rb') as f:
            genome = pickle.load(f)

        full_readcounts = {chrom: [0 for i in range(num_regions[chrom])] for chrom in chrom_names}
        for chrom in chrom_names:
            for allele in [0, 1]:
                if len(genome[chrom][allele]) > 0:
                    for homolog in genome[chrom][allele]:
                        while 'X' in homolog:
                            homolog.remove('X')
                        
                        num_windows = int(np.ceil(len(homolog) / regions_per_window))
                        if uniform_coverage:
                            cur_readcounts = [round(poisson.rvs((window_size*coverage) / (2*read_len))) for i in range(num_windows)]
                        else:
                            cur_readcounts = [round(x) for x  in draw_readcounts(num_windows, window_size, interval, Aa, Bb, coverage, read_len)]

                        windows = [list(w) for w in iter_by_chunk(homolog, regions_per_window)]
                        while None in windows[-1]:
                            windows[-1].remove(None)

                        for i,window in enumerate(windows):
                            rc = cur_readcounts[i]
                            cur_window_len = len(window)
                            if cur_window_len != regions_per_window:
                                rc = round(rc * (cur_window_len/regions_per_window))
                            
                            counts = dict(Counter(window))
                            for r,c in counts.items():
                                r_ratio = c / cur_window_len
                                full_readcounts[chrom][r] += (r_ratio * rc)
        

        binned_readcounts = {chrom: [] for chrom in chrom_names}
        for chrom in chrom_names:
            for i in range(len(bins[chrom])-1):
                regions = list(range(bins[chrom][i], bins[chrom][i+1]))
                binned_readcounts[chrom].append(round(sum([full_readcounts[chrom][j] for j in regions])))
        final_readcounts[cell] = binned_readcounts
        os.remove(prefix + '.pkl')

    with open(os.path.join(out_path, 'readcounts.tsv'), 'w+') as f:
        headers = ['CELL', 'chrom', 'start', 'end', 'readcount']
        f.write('\t'.join(headers) + '\n')

        # Writes to file in order of chroms, bins, cell
        for chrom in chrom_names:
            num_bins = len(bins[chrom]) - 1
            for i in range(num_bins):
                bin_start, bin_end = str(bins[chrom][i]*region_length), str(bins[chrom][i+1]*region_length)
                for cell in leaves:
                    line = [cell.name, chrom, bin_start, bin_end, str(final_readcounts[cell][chrom][i])]
                    f.write('\t'.join(line) + '\n')




    