import subprocess
import os

#call = subprocess.run(['dwgsim', '-H', '-N', '10', '-o', '1', '/Users/samsonweiner/Desktop/Hippo/data/chr20_subset1.fa', 'test1'])

from pyfaidx import Fasta
from sequence import *
from reads import *
from tree import Node

genome = {}
genome['chr20'] = [[], []]
genome['chr20'][0] = [i for i in range(1000)]
genome['chr20'][1] = [i for i in range(1000)]

genome['chr20'][0] = genome['chr20'][0][:300] + 2*genome['chr20'][0][200:300] + genome['chr20'][0][300:]
genome['chr20'][1] = genome['chr20'][1][:600] + genome['chr20'][1][750:]


ref, ref_chrom_lens = read_fasta('/Users/samsonweiner/Desktop/Hippo/data/chr20_subset1.fa')
#ref, ref_chrom_lens = read_fasta('/Users/samsonweiner/Desktop/Hippo/reference/hg38-chroms/chr20.fa')
#print(ref_chrom_lens['chr20']/1000)

[Aa, Bb] = get_alpha_beta(0.5, 0.28)

c = Node(name='cell1')
c.genome = genome

gen_reads_cell(c, ref, ['chr20'], 1000, 200000, 3, Aa, Bb, 0.2, 35)



#allele = 0
#chrom = 'chr20'
#cell_ref, chrom_lens = build_cell_ref(c.genome, ref, 1000, allele, 'cell1')

#window_size = 200000
#num_windows = round(chrom_lens[chrom] / window_size)
#interval = 3
#coverage = 0.2
#readlen = 35

#readcounts = draw_readcounts(num_windows, window_size, interval, Aa, Bb, coverage, readlen)

#w = 1
#start = w*window_size + 1
#end = (w+1)*window_size
#region_fa_path = cell_ref[:-3] + '_region' + str(start) + '.fa'

#f = open(region_fa_path, 'w+')
#call = subprocess.run(['samtools', 'faidx', cell_ref, chrom + ':' + str(start) + '-' + str(end)], stdout=f)
#f.close()

#call = subprocess.run(['dwgsim', '-H', '-o', '1', '-N', str(readcounts[w]), '-1', str(readlen), '-2', str(readlen), region_fa_path, region_fa_path[:-3]])

#os.rename(region_fa_path[:-3] + '.bwa.read1.fastq.gz', cell_ref[:-3] + '.read1.fastq.gz')
#os.rename(region_fa_path[:-3] + '.bwa.read2.fastq.gz', cell_ref[:-3] + '.read2.fastq.gz')

#with open(cell_ref[:-3] + '.read1.fastq.gz', 'a') as f1, open(cell_ref[:-3] + '.read2.fastq.gz', 'a') as f2:
#    call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read1.fastq.gz'], stdout=f1)
#    call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read2.fastq.gz'], stdout=f2)

#f = open(cell_ref[:-3] + '.read2.fastq.gz', 'a')
#call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read2.fastq.gz'], stdout=f)
#f.close()

#call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read1.fastq.gz', '>>', cell_ref[:-3] + '.read1.fastq.gz'])
#call = subprocess.run(['cat', region_fa_path[:-3] + '.bwa.read2.fastq.gz', '>>', cell_ref[:-3] + '.read2.fastq.gz'])
#call = subprocess.run(['rm', region_fa_path[:-3] + '.bwa.read1.fastq.gz', region_fa_path[:-3] + '.bwa.read2.fastq.gz'])


#call = subprocess.run(['rm', region_fa_path, region_fa_path[:-3] + '.mutations.txt', region_fa_path[:-3] + '.mutations.vcf'])