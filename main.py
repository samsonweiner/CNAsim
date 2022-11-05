import argparse
import os
import time

from tree import *
from sim_genomes import *
from noise import *
from format_profiles import *
from utilities import *

def parse_args():
    #Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--out-path', type=str, default='/', help='Path to output directory.')
    parser.add_argument('-t', '--tree-type', type=int, default=0, help='0: ms, 1: random, 2: from file (replace ms-path)')
    parser.add_argument('-s', '--growth-rate', type=float, default=15.1403, help='Exponential growth rate for ms.')
    parser.add_argument('-f', '--ms-path', type=str, default='../../msdir/ms', help='Path to ms binary.')
    parser.add_argument('-n', '--num-cells', type=int, default=250, help='Number of observed cells in sample.')
    parser.add_argument('-p', '--placement-type', type=int, default=0, help='0: draw from a Poisson with fixed mean, 1: draw from a Poisson with mean prop to edge length, 2: fixed per edge')
    parser.add_argument('-c', '--placement-param', type=float, default=2, help='Parameter for placement choice.')
    parser.add_argument('-e', '--min-cn-length', type=int, default=1000, help='Minimum copy number event length in bp.')
    parser.add_argument('-l', '--cn-length-mean', type=int, default=5000000, help='Mean copy number event length in bp.')
    parser.add_argument('-m', '--cn-copy-param', type=float, default=0.5, help='Parameter in the geometric to select number of copies.')
    parser.add_argument('-g', '--cn-event-rate', type=float, default=0.5, help='Probability an event is an amplification. Deletion rate is 1 - amp rate.')
    parser.add_argument('-r', '--root-event-mult', type=int, default=10, help='Multiplier for the number of events along edge into founder cell.')
    parser.add_argument('-w', '--whole-genome-dup', action='store_true', help='Include WGD.')
    parser.add_argument('-v', '--whole-chrom-event', action='store_true', help='Include whole chromosomal alterations.')
    parser.add_argument('-i', '--whole-chrom-rate', type=float, default=0.2, help='Probability that a chromosome is affected by a whole chrom event.')
    parser.add_argument('-u', '--whole-chrom-type', type=float, default=0.5, help='Probability that a whole chrom event is an amplification.')
    parser.add_argument('-q', '--whole-chrom-copy', type=float, default=0.8, help='Parameter in the geometric to select the number of copies for whole chrom events.')
    parser.add_argument('-x', '--num-chromosomes', type=int, default=22, help='Number of chromosomes.')
    parser.add_argument('-y', '--chrom-length', type=int, default=100000000, help='Length of chromosomes in bp.')
    parser.add_argument('-z', '--use-hg38-lengths', action='store_true', help='Use hg38 chromosome lengths.')
    parser.add_argument('-b', '--bin-length', type=int, default=5000000, help='Resolution of copy number profiles in bp.')
    parser.add_argument('-a', '--error-type', type=int, default=0, help='0: none, 1: simple, 2: normal, 3: boundary, 4: sequence, 5: mixed')
    parser.add_argument('-a1', '--error-rate-1', type=float, default=0, help='Primary parameter for introducing error.')
    parser.add_argument('-a2', '--error-rate-2', type=float, default=0, help='Secondary parameter for introducing error.')
    parser.add_argument('-d', '--summary', action='store_true', help='Summarize simulation statistics.')
    parser.add_argument('-k', '--param-file', type=str, default=None, help='Optional parameter file.')
    arguments = parser.parse_args()

    return handle_args(arguments)

def main(args):
    #Output dir
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
    if args['tree_type'] == 0:
        tree_str = call_ms(args['ms_path'], args['num_cells'], args['out_path'], args['growth_rate'])
        tree = Tree(newick=tree_str)
        tree.set_node_names()
        tree.save(os.path.join(args['out_path'], 'tree.nwk'), format=1)

    elif args['tree_type'] == 1:
        tree = gen_random_topology(args['num_cells'])
        tree.set_node_names()
        tree.save(os.path.join(args['out_path'], 'tree.nwk'))
    
    else:
        f1 = open(args['ms_path'])
        tree_str = f1.readline().strip()
        f1.close()
        tree = Tree(newick=tree_str)
        tree.set_node_names()
        tree.save(os.path.join(args['out_path'], 'tree.nwk'), format=1)

    #Simulate evolution
    print('Generating genomes, events, and profiles...')
    chrom_names = ['chr' + str(i+1) for i in range(args['num_chroms'])]
    normal_diploid_genome = init_diploid_genome(args['min_cn_length'], chrom_names, args['chrom_length'], args['use_hg38_lengths'])
    bins = gen_profiles(args, tree, chrom_names, normal_diploid_genome)

    #Format and output profiles
    print('Formating and saving profiles')
    if args['error_type'] != 0:
        save_CN_profiles(tree, chrom_names, bins, args['min_cn_length'], os.path.join(args['out_path'], 'clean_profiles.tsv'))
        choose_noise_model(args['error_type'], tree, chrom_names, args['error_rate_1'], args['error_rate_2'])
    save_CN_profiles(tree, chrom_names, bins, args['min_cn_length'], os.path.join(args['out_path'], 'profiles.tsv'))

    #Log summary stats and execution time
    total_time = round(time.time() - start)
    log.write('\nTime: ' + str(total_time) + '\n')
    log.close()
    if args['summary']:
        summary(tree, os.path.join(args['out_path'], 'summary.txt'), args['num_chroms'], args['WGD'], args['min_cn_length'])


if __name__ == '__main__':
    args = parse_args()
    main(args)