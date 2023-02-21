import argparse
import os
import time
import sys

from tree import *
from sim_genomes import *
from sequence import *
from reads import *
from noise import *
from format_profiles import *
from utilities import *

def parse_args():
    #Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', type=int, required=True, help='Main simulator mode for generating data. 0: CNP data, 1: read data, 2: both')
    parser.add_argument('-o', '--out-path', type=str, default='/', help='Path to output directory.')
    parser.add_argument('-t', '--tree-type', type=int, default=0, help='0: ms, 1: random, 2: from file (use -T to specify file path).')
    parser.add_argument('-T', '--tree-path', type=str, default='', help='Path to input tree.')
    parser.add_argument('-g', '--growth-rate', type=float, default=15.1403, help='Exponential growth rate for ms.')
    parser.add_argument('-f', '--ms-path', type=str, default='ms', help='Path to ms binary.')
    parser.add_argument('-n', '--num-cells', type=int, default=250, help='Number of observed cells in sample.')
    parser.add_argument('-n1', '--normal-fraction', type=float, default=0, help='Proportion of cells that are normal.')
    parser.add_argument('-n2', '--pseudonormal-fraction', type=float, default=0, help='Proportion of cells that are pseudonormal cells.')
    parser.add_argument('-c', '--num-clones', type=int, default=0, help='Number of ancestral nodes to select as clonal founders.')
    parser.add_argument('-p1', '--placement-type', type=int, default=0, help='0: draw from a Poisson with fixed mean, 1: draw from a Poisson with mean prop to edge length, 2: fixed per edge')
    parser.add_argument('-p2', '--placement-param', type=float, default=2, help='Parameter for placement choice.')
    parser.add_argument('-k', '--min-cn-length', type=int, default=1000, help='Minimum copy number event length in bp.')
    parser.add_argument('-l', '--cn-length-mean', type=int, default=5000000, help='Mean copy number event length in bp.')
    parser.add_argument('-a', '--cn-copy-param', type=float, default=0.5, help='Parameter in the geometric to select number of copies.')
    parser.add_argument('-s', '--cn-event-rate', type=float, default=0.5, help='Probability an event is an amplification. Deletion rate is 1 - amp rate.')
    parser.add_argument('-j', '--founder-event-mult', type=int, default=10, help='Multiplier for the number of events along edge into founder cell.')

    parser.add_argument('-w', '--WGD', action='store_true', help='Include WGD.')
    parser.add_argument('-v', '--chrom-level-event', action='store_true', help='Include whole chromosomal alterations.')
    parser.add_argument('-q', '--chrom-arm-rate', type=float, default=0.75, help='Probability that a chromosomal event is a chromosome-arm event.')
    parser.add_argument('-i1', '--chrom-rate-founder', type=int, default=2, help='Parameter in poisson for number of chrom-level events along the edges into and out of the founder cell.')
    parser.add_argument('-i2', '--chrom-rate-clone', type=int, default=1.5, help='Parameter in poisson for number of chrom-level events for clonal nodes.')
    parser.add_argument('-u', '--chrom-event-type', type=float, default=0.5, help='Probability that a whole chrom event is a duplication.')

    parser.add_argument('-N', '--num-chromosomes', type=int, default=22, help='Number of chromosomes.')
    parser.add_argument('-L', '--chrom-length', type=int, default=100000000, help='Length of chromosomes in bp.')
    parser.add_argument('-A', '--chrom-arm-ratio', type=float, default=0.5, help='If fixed size chromosomes, ratio of length within the p-arm.')
    parser.add_argument('-B', '--bin-length', type=int, default=5000000, help='Resolution of copy number profiles in bp.')
    parser.add_argument('-E1', '--error-rate-1', type=float, default=0, help='Error rate for the boundary model.')
    parser.add_argument('-E2', '--error-rate-2', type=float, default=0, help='Error rate for the jitter model.')
    parser.add_argument('-U', '--use-hg38-static', action='store_true', help='Use hg38 chromosome information.')

    parser.add_argument('-r', '--reference', type=str, default='', help='Path to input reference genome in fasta format.')
    parser.add_argument('-M', '--use-uniform-coverage', action='store_true', help='Use uniform coverage across the genome.')
    parser.add_argument('-X', '--lorenz-x', type=float, default=0.5, help='x-coordinate for point on lorenz curve')
    parser.add_argument('-Y', '--lorenz-y', type=float, default=0.4, help='y-coordinate for point on lorenz curve')
    parser.add_argument('-W', '--window-size', type=int, default=1000000, help='Number of base pairs to generate reads for in each iteration.')
    parser.add_argument('-I', '--interval', type=int, default=3, help='Initializes a point in the coverage distribution every interval number of windows.')
    parser.add_argument('-C', '--coverage', type=float, default=0.1, help='Average sequencing coverage across the genome.')
    parser.add_argument('-R', '--read-length', type=int, default=100, help='Paired-end short read length.')
    parser.add_argument('-P', '--processors', type=int, default=1, help='Number of processes to use for generating reads in parallel.')

    parser.add_argument('-d', '--summary', action='store_true', help='Summarize simulation statistics.')
    parser.add_argument('-F', '--param-file', type=str, default=None, help='Optional parameter file.')
    arguments = parser.parse_args()

    return handle_args(arguments)

def main(args):
    # Initialize output directory
    print('Output directory', os.path.abspath(args['out_path']))
    if not os.path.isdir(args['out_path']):
        os.makedirs(args['out_path'])

    #Log parameters
    start = time.time()
    log = open(os.path.join(args['out_path'], 'log.txt'), 'w+')
    for k,v in args.items():
        log.write(k + ': ' + str(v) + '\n')

    #Create tree structure
    print('Preparing ground truth tree...')
    founder_events = args['placement_param'] * args['founder_event_mult']
    tree = make_tumor_tree(args['tree_type'], args['num_cells'], args['normal_fraction'], args['pseudonormal_fraction'], founder_events, args['out_path'], args['ms_path'], args['growth_rate'], args['tree_path'])
    if args['num_clones'] > 0:
        select_clones(tree, args['num_clones'])
    tree.set_node_names()
    tree.save(os.path.join(args['out_path'], 'tree.nwk'), format=3)

    record_cell_types(tree, os.path.join(args['out_path'], 'cell_types.tsv'))

    #if args['placement_type'] == 1:
    #    scale_edge_lengths(tree, args['placement_param'])

    #Simulate evolution
    print('Generating genomes, events, and profiles...')
    if args['use_hg38_static']:
        file_path = os.path.realpath(__file__)
        sim_dir_path, filename = os.path.split(file_path)
        chrom_lens, arm_ratios = hg38_chrom_lengths_from_cytoband(os.path.join(sim_dir_path, 'resources/cytoBand.txt'), include_allosomes=False, include_arms=True)
    else:
        chrom_lens = {}
        for i in range(1, args['num_chromosomes']+1):
            chrom_lens['chr' + str(i)] = args['chrom_length']
    
    '''
    # Sequence mode
    if args['mode'] == 1:
        ref, ref_chrom_lens = read_fasta(args['reference'])
        ref_chrom_lens.pop('chrX', None)
        ref_chrom_lens.pop('chrY', None)
        chrom_names = list(ref_chrom_lens.keys())

        if args['use_hg38_static']:
            normal_diploid_genome, num_regions = init_diploid_genome(args['min_cn_length'], chrom_names, ref_chrom_lens, arm_ratios)
        else:
            normal_diploid_genome, num_regions = init_diploid_genome(args['min_cn_length'], chrom_names, ref_chrom_lens, args['chrom_arm_ratio'])
        tree.root.genome = normal_diploid_genome

        evolve_tree(tree.root, args, chrom_names, num_regions)

        print('Generating reads')
        gen_reads(ref, num_regions, chrom_names, tree, args['use_uniform_coverage'], args['lorenz_x'], args['lorenz_y'], args['min_cn_length'], args['interval'], args['window_size'], args['coverage'], args['read_length'], args['out_path'], args['processors'])

    # CNP mode
    if args['mode'] == 0:
        chrom_names = ['chr' + str(i+1) for i in range(args['num_chromosomes'])]
        if args['use_hg38_static']:
            normal_diploid_genome, num_regions = init_diploid_genome(args['min_cn_length'], chrom_names, chrom_lens, arm_ratios)
        else:
            normal_diploid_genome, num_regions = init_diploid_genome(args['min_cn_length'], chrom_names, args['chrom_length'], args['chrom_arm_ratio'])
        tree.root.genome = normal_diploid_genome

        regions_per_bin = np.floor(args['bin_length']/args['min_cn_length'])
        bins = {}
        for chrom in chrom_names:
            bins[chrom] = [k for k in range(num_regions[chrom]) if k % regions_per_bin == 0]
            if (num_regions[chrom] - 1) - bins[chrom][-1] < regions_per_bin/2:
                bins[chrom][-1] = num_regions[chrom]
            else:
                bins[chrom].append(num_regions[chrom])
            
        evolve_tree(tree.root, args, chrom_names, num_regions, bins=bins)

        print('Formating and saving profiles')
        if args['error_rate_1'] != 0 or args['error_rate_2'] != 0:
            #save_CN_profiles(tree, chrom_names, bins, args['min_cn_length'], os.path.join(args['out_path'], 'clean_profiles.tsv'))
            add_noise_mixed(tree, chrom_names, args['error_rate_1'], args['error_rate_2'])
        save_CN_profiles(tree, chrom_names, bins, args['min_cn_length'], os.path.join(args['out_path'], 'profiles.tsv'))
    '''

    ##########################
    if os.path.isfile(args['reference']):
        ref, chrom_lens = read_fasta(args['reference'])
        chrom_lens.pop('chrX', None)
        chrom_lens.pop('chrY', None)
        chrom_names = list(chrom_lens.keys())
    else:
        chrom_names = ['chr' + str(i+1) for i in range(args['num_chromosomes'])]

    if args['use_hg38_static']:
        normal_diploid_genome, num_regions = init_diploid_genome(args['min_cn_length'], chrom_names, chrom_lens, arm_ratios)
    else:
        normal_diploid_genome, num_regions = init_diploid_genome(args['min_cn_length'], chrom_names, chrom_lens, args['chrom_arm_ratio'])
    tree.root.genome = normal_diploid_genome

    if args['mode'] == 0 or args['mode'] == 2:
        regions_per_bin = np.floor(args['bin_length']/args['min_cn_length'])
        bins = {}
        for chrom in chrom_names:
            bins[chrom] = [k for k in range(num_regions[chrom]) if k % regions_per_bin == 0]
            if (num_regions[chrom] - 1) - bins[chrom][-1] < regions_per_bin/2:
                bins[chrom][-1] = num_regions[chrom]
            else:
                bins[chrom].append(num_regions[chrom])
    else:
        bins = None

    evolve_tree(tree.root, args, chrom_names, num_regions, bins=bins)

    if args['mode'] == 0 or args['mode'] == 2:
        print('Formating and saving profiles')
        if args['error_rate_1'] != 0 or args['error_rate_2'] != 0:
            #save_CN_profiles(tree, chrom_names, bins, args['min_cn_length'], os.path.join(args['out_path'], 'clean_profiles.tsv'))
            add_noise_mixed(tree, chrom_names, args['error_rate_1'], args['error_rate_2'])
        save_CN_profiles(tree, chrom_names, bins, args['min_cn_length'], os.path.join(args['out_path'], 'profiles.tsv'))
    
    if args['mode'] == 1 or args['mode'] == 2:
        gen_reads(ref, num_regions, chrom_names, tree, args['use_uniform_coverage'], args['lorenz_x'], args['lorenz_y'], args['min_cn_length'], args['interval'], args['window_size'], args['coverage'] / 2, args['read_length'], args['out_path'], args['processors'])
    ###########################

    #Log summary stats and execution time
    total_time = round(time.time() - start)
    log.write('\nTime: ' + str(total_time) + '\n')
    log.close()
    #if args['summary']:
        #summary(tree, os.path.join(args['out_path'], 'summary.txt'), args['num_chroms'], args['WGD'], args['min_cn_length'])
    #if args['mode'] == 0 or args['mode'] == 1:
    record_events(tree, os.path.join(args['out_path'], 'events.tsv'))

if __name__ == '__main__':
    #profiler = cProfile.Profile()
    #profiler.enable()
    args = parse_args()
    main(args)
    #profiler.disable()
    #stats = pstats.Stats(profiler).sort_stats('cumtime')
    #stats.dump_stats(os.path.join(args['out_path'], 'stats.txt'))