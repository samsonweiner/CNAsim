import os
import subprocess

import numpy as np
from scipy.stats import beta
from scipy.optimize import newton_krylov
from scipy.optimize.nonlin import NoConvergence

from sequence import *

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
    new_points = [max(min(i[1], 1), 0) for i in new_points]
    return new_points

# Draw readcounts for each bin
def draw_readcounts(num_windows, window_size, interval, Aa, Bb, coverage, readlen):
    avg_read = (coverage*window_size) / readlen
    cov_scales = gen_coverage(num_windows, interval, Aa, Bb)
    exp_counts = [2*avg_read*x for x in cov_scales]
    readcounts = [np.random.poisson(x) for x in exp_counts]
    return readcounts


# Generate reads across the genome for a given cell
def gen_reads_cell(cell, ref, chrom_names, min_cn_len, window_size, interval, Aa, Bb, coverage, readlen):
    for allele in [0, 1]:
        cell_ref, chrom_lens = build_cell_ref(cell.genome, ref, min_cn_len, allele, cell.name)
        init = False
        for chrom in chrom_names:
            num_windows = round(chrom_lens[chrom] / window_size)
            readcounts = draw_readcounts(num_windows, window_size, interval, Aa, Bb, coverage, readlen)
            for w in range(num_windows):
                start = w * window_size + 1
                if w == num_windows - 1:
                    end = chrom_lens[chrom]
                else:
                    end = (w+1) * window_size
                
                region_fa_path = cell_ref[:-3] + '_region' + str(start) + '.fa'
                with open(region_fa_path, 'w+') as f:
                    call = subprocess.run(['samtools', 'faidx', cell_ref, chrom + ':' + str(start) + '-' + str(end)], stdout=f)

                proc = subprocess.run(['dwgsim', '-H', '-o', '1', '-N', str(readcounts[w]), '-1', str(readlen), '-2', str(readlen), region_fa_path, region_fa_path[:-3]])
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
                os.remove(region_fa_path)
                os.remove(region_fa_path[:-3] + '.mutations.txt')
                os.remove(region_fa_path[:-3] + '.mutations.vcf') 
        os.remove(cell_ref)
        os.remove(cell_ref + '.fai')

# Generate reads for all cells
def gen_reads(ref_regions, tree, x0, y0, interval, window, coverage, readlen):
    [Aa, Bb] = get_alpha_beta(x0, y0)
