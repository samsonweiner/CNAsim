#!/usr/bin/env python3

# Copyright (C) 2023 Samson Weiner (samson.weiner@uconn.edu) and
# Mukul S. Bansal (mukul.bansal@uconn.edu).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.


import argparse
import os
import time

from .tree import *
from .sim_genomes import init_diploid_genome, evolve_tree, prepare_ancestral_profiles
from .sequence import read_fasta
from .reads import gen_reads, gen_readcounts
from .noise import add_noise_mixed
from .format_profiles import save_CN_profiles_leaves, save_CN_profiles_ancestors
from .utilities import *


def parse_args():
    #Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mode', type=int, default=None, help='Main simulator mode for generating data. 0: CNPs only, 1: readcounts & CNPs, 2: sequencing reads & CNPs')
    parser.add_argument('-o', '--out-path', type=str, default='CNAsim_output/', help='Path to output directory.')
    parser.add_argument('-t', '--tree-type', type=int, default=0, help='0: coalescence, 1: random, 2: from file (use -T to specify file path).')
    parser.add_argument('-T', '--tree-path', type=str, default=None, help='Path to input tree.')
    parser.add_argument('-g', '--growth-rate', type=float, default=0.003785, help='Exponential growth rate for standard coalescent.')
    parser.add_argument('-s', '--num-sweep', type=int, default=0, help='Number of selective sweeps in coalescent.')
    parser.add_argument('-s1', '--selection-strength', type=float, default=0.01, help='Parameter controlling the strength of selection during the sweeps.')
    parser.add_argument('-n', '--num-cells', type=int, default=250, help='Number of observed cells in sample.')
    parser.add_argument('-n1', '--normal-fraction', type=float, default=0, help='Proportion of cells that are normal.')
    parser.add_argument('-n2', '--pseudonormal-fraction', type=float, default=0, help='Proportion of cells that are pseudonormal cells.')
    parser.add_argument('-c', '--num-clones', type=int, default=0, help='Number of ancestral nodes to select as clonal founders.')
    parser.add_argument('-c1', '--clone-criteria', type=int, default=0, help='Criteria to choose clonal ancesters. 0: proportional to number of leaves in subtree. 1: proportional to edge length')
    parser.add_argument('-c2', '--clone-mu', type=float, default=None, help='Mean number of leaves in a subclone. Must select 0 for clone-criteria.')
    parser.add_argument('-c3', '--clone-sd', type=float, default=None, help='SD in the number of leaves in a subclone. Must select 0 for clone-criteria.')
    parser.add_argument('-p', '--placement-type', type=int, default=0, help='Number of CNAs per edge. 0: draw from a Poisson with fixed mean, 1: draw from a Poisson with mean prop to edge length, 2: fixed per edge')
    parser.add_argument('-p1', '--placement-param', type=float, default=2, help='Parameter for placement choice.')
    parser.add_argument('-k', '--region-length', type=int, default=1000, help='Region length in bp. Essentially controls the resolution of the simulated genome.')
    parser.add_argument('-l', '--cn-length-mean', type=int, default=5000000, help='Mean copy number event length in bp.')
    parser.add_argument('-l1', '--min-cn-length', type=int, default=1000, help='Minimum copy number event length in bp. Should be at minimum the region length and less than the mean.')
    parser.add_argument('-a', '--cn-copy-param', type=float, default=0.5, help='Parameter in the geometric to select number of copies.')
    parser.add_argument('-b', '--cn-event-rate', type=float, default=0.5, help='Probability an event is an amplification. Deletion rate is 1 - amp rate.')
    parser.add_argument('-j', '--founder-event-mult', type=int, default=10, help='Multiplier for the number of events along edge into founder cell.')
    parser.add_argument('-w', '--WGD', action='store_true', help='Include WGD.')
    parser.add_argument('-v', '--chrom-level-event', action='store_true', help='Include chromosomal alterations.')
    parser.add_argument('-q', '--chrom-arm-rate', type=float, default=0.75, help='Probability that a chromosomal event is a chromosome-arm event.')
    parser.add_argument('-i', '--chrom-rate-founder', type=float, default=2, help='Parameter in poisson for number of chromosome-level events along the edge into the founder cell.')
    parser.add_argument('-i1', '--chrom-rate-super-clone', type=float, default=1, help='Parameter in poisson for number of chromosome-level events along the edges out of the founder cell.')
    parser.add_argument('-i2', '--chrom-rate-clone', type=float, default=1, help='Parameter in poisson for number of chrom-level events for clonal nodes.')
    parser.add_argument('-u', '--chrom-event-type', type=float, default=0.5, help='Probability that a chromosomal event is a duplication.')
    parser.add_argument('-N', '--num-chromosomes', type=int, default=22, help='Number of chromosomes if run in mode 0.')
    parser.add_argument('-L', '--chrom-length', type=int, default=100000000, help='Length of chromosomes in bp if not using hg38 static.')
    parser.add_argument('-A', '--chrom-arm-ratio', type=float, default=0.5, help='If not using hg38 static, ratio of length within the p-arm.')
    parser.add_argument('-B', '--bin-length', type=int, default=1000000, help='Resolution of copy number profiles in bp.')
    parser.add_argument('-E1', '--error-rate-1', type=float, default=0, help='Error rate for the boundary model.')
    parser.add_argument('-E2', '--error-rate-2', type=float, default=0, help='Error rate for the jitter model.')
    #parser.add_argument('-O', '--output-clean-CNP', action='store_true', help='Output the clean CNPs in addition to the noisy ones.')
    parser.add_argument('-U', '--use-hg38-static', action='store_true', help='Use hg38 chromosome information. Excludes sex chromosomes chrX and chrY.')
    #parser.add_argument('-A', '--include-autosomes', action='store_true', help='Include autosomes in genome. Must either use hg38 static with number of chromosomes > 22, or pass a reference genome with autosomes.')
    parser.add_argument('-r1', '--reference', type=str, default=None, help='Path to input reference genome as the primary haplotype in fasta format. Will be duplicated as both haplotypes if an alternate is not provided.')
    parser.add_argument('-r2', '--alt-reference', type=str, default=None, help='Path to an alternate reference genome to be used as a secondary haplotype, also in fasta format.')
    #parser.add_argument('-D', '--chrom-names-path', type=str, default=None, help='Path to file containing chromosomes names to use in reference file(s). By default, will use all chromosomes present in reference except chrX and chrY.')
    parser.add_argument('-M', '--use-uniform-coverage', action='store_true', help='Use uniform coverage across the genome.')
    parser.add_argument('-X', '--lorenz-x', type=float, default=0.5, help='x-coordinate of point on lorenz curve if using non-uniform coverage.')
    parser.add_argument('-Y', '--lorenz-y', type=float, default=0.4, help='y-coordinate of point on lorenz curve if using non-uniform coverage.')
    parser.add_argument('-W', '--window-size', type=int, default=1000000, help='Number of base pairs to generate reads for in each iteration.')
    parser.add_argument('-I', '--interval', type=int, default=3, help='Initializes a point in the coverage distribution every interval number of windows.')
    parser.add_argument('-C', '--coverage', type=float, default=0.1, help='Average sequencing coverage across the genome.')
    parser.add_argument('-R', '--read-length', type=int, default=150, help='Paired-end short read length.')
    parser.add_argument('-S', '--seq-error', type=float, default=0.02, help='Per base error rate for generating sequence data.')
    parser.add_argument('-P', '--processors', type=int, default=1, help='Number of processes to use for generating reads in parallel.')
    parser.add_argument('-d', '--disable-info', action='store_true', help='Do not output simulation log, cell types, or ground truth events.')
    parser.add_argument('-F', '--param-file', type=str, default=None, help='Path to parameter file to specify parameters instead of the command line. Must conform to the sample format.')
    arguments = parser.parse_args()

    return handle_args(arguments)

def main():
    args = parse_args()

    if args['mode'] not in [0, 1, 2]:
        raise ModeError

    ## Initialize output directory
    print('Output directory:', os.path.abspath(args['out_path']))
    if not os.path.isdir(args['out_path']):
        os.makedirs(args['out_path'])

    ## Log parameters
    if not args['disable_info']:
        start = time.time()
        log = open(os.path.join(args['out_path'], 'log.txt'), 'w+')
        for k,v in args.items():
            log.write(k + ': ' + str(v) + '\n')

    ## Create tree structure
    print('Preparing ground truth tree...')
    founder_events = int(args['placement_param'] * args['founder_event_mult'])
    tree = make_tumor_tree(args['tree_type'], args['num_cells'], args['normal_fraction'], args['pseudonormal_fraction'], founder_events, args['out_path'], args['growth_rate'], args['tree_path'], args['num_sweep'], args['selection_strength'])
    clone_founders = []
    if args['num_clones'] > 0:
        clone_founders = select_clones(tree, args['num_clones'], args['clone_criteria'], args['clone_mu'], args['clone_sd'])
        args['chrom_level_event'] = True
    tree.set_node_names()
    tree.save(os.path.join(args['out_path'], 'tree.nwk'), format=3)

    if not args['disable_info']:
        clone_founders = record_cell_types(tree, os.path.join(args['out_path'], 'cell_types.tsv'), args['chrom_level_event'], args['chrom_rate_super_clone'], clone_founders=clone_founders)

    if args['placement_type'] == 1:
        scale_edge_lengths(tree, args['placement_param'])

    ## Initialize genome
    print('Initializing reference genome...')
    if args['reference']:
        if not os.path.isfile(args['reference']):
            raise InputError("Cannot access input " + str(args['reference']))
        ref1, chrom_lens = read_fasta(args['reference'])
        if args['use_hg38_static']:
            all_chroms = list(chrom_lens.keys())
            chrom_names = [f'chr{i}' for i in range(1, 23)]
            for c in all_chroms:
                if c not in chrom_names:
                    chrom_lens.pop(c, None)
        chrom_names = list(chrom_lens.keys())
        print(f'Number of chromosomes in reference: {len(chrom_names)}')
        if not check_chrom_lengths(chrom_lens, args['window_size']):
            if args['mode'] == 2:
                raise ChromSizeError()

        if args['alt_reference']:
            if not os.path.isfile(args['alt_reference']):
                raise InputError("Cannot access input " + str(args['alt_reference']))
            ref2, alt_chrom_lens = read_fasta(args['alt_reference'])
            if args['use_hg38_static']:
                all_chroms = list(chrom_lens.keys())
                alt_chrom_names = [f'chr{i}' for i in range(1, 23)]
                for c in all_chroms:
                    if c not in alt_chrom_names:
                        alt_chrom_lens.pop(c, None)
            if set(chrom_names) != set(alt_chrom_lens.keys()):
                raise ChromNameError()
        else:
            ref2 = ref1
    else:
        if args['mode'] == 2:
            raise InputError("Cannot access input " + str(args['reference']))
        chrom_names = ['chr' + str(i+1) for i in range(args['num_chromosomes'])]
        chrom_lens = dict(zip(chrom_names, [args['chrom_length'] for i in range(args['num_chromosomes'])]))

    if args['use_hg38_static']:
        file_path = os.path.realpath(__file__)
        sim_dir_path, filename = os.path.split(file_path)
        if not args['reference']:
            chrom_lens, arm_ratios = hg38_chrom_lengths_from_cytoband(include_allosomes=False, include_arms=True)
            if args['num_chromosomes'] > 22:
                raise ChromNumError()
        else:
            temp, arm_ratios = hg38_chrom_lengths_from_cytoband(include_allosomes=False, include_arms=True)
            for chrom in chrom_names:
                if chrom not in arm_ratios:
                    raise ChromNameError()
    else:
        arm_ratios = args['chrom_arm_ratio']
    
    normal_diploid_genome, num_regions = init_diploid_genome(args['region_length'], chrom_names, chrom_lens, arm_ratios)
    tree.root.genome = normal_diploid_genome

    ## Evolve genomes along tree
    print('Generating genomes, events, and profiles...')
    #if args['mode'] == 0 or args['mode'] == 2:
    regions_per_bin = np.floor(args['bin_length']/args['region_length'])
    bins = {}
    for chrom in chrom_names:
        bins[chrom] = [k for k in range(num_regions[chrom]) if k % regions_per_bin == 0]
        if (num_regions[chrom] - 1) - bins[chrom][-1] < regions_per_bin/2:
            bins[chrom][-1] = num_regions[chrom]
        else:
            bins[chrom].append(num_regions[chrom])

    evolve_tree(tree.root, args, chrom_names, num_regions, bins=bins)
    prepare_ancestral_profiles(tree, args, chrom_names, num_regions, bins=bins)

    ## Generate data
    #if args['mode'] == 0 or args['mode'] == 2:
    print('Formating and saving profiles')
    if args['error_rate_1'] != 0 or args['error_rate_2'] != 0:
        save_CN_profiles_leaves(tree, chrom_names, bins, args['region_length'], os.path.join(args['out_path'], 'clean_profiles.tsv'))
        add_noise_mixed(tree, chrom_names, args['error_rate_1'], args['error_rate_2'])
    save_CN_profiles_leaves(tree, chrom_names, bins, args['region_length'], os.path.join(args['out_path'], 'profiles.tsv'))
    save_CN_profiles_ancestors(tree, chrom_names, bins, args['region_length'], os.path.join(args['out_path'], 'ancestral_profiles.tsv'))

    if args['mode'] == 1:
        print('Generating read counts...')
        gen_readcounts(tree, chrom_names, bins, num_regions, args['region_length'], args['use_uniform_coverage'], args['lorenz_x'], args['lorenz_y'], args['interval'], args['window_size'], args['coverage'] / 2, args['read_length'], args['out_path'])
    
    if args['mode'] == 2:
        print('Generating synthetic sequencing reads...')
        gen_reads(ref1, ref2, num_regions, chrom_names, tree, args['use_uniform_coverage'], args['lorenz_x'], args['lorenz_y'], args['region_length'], args['interval'], args['window_size'], args['coverage'] / 2, args['read_length'], args['seq_error'], args['out_path'], args['processors'])

    ## Logging information
    if not args['disable_info']:
        total_time = round(time.time() - start)
        log.write('\nTime: ' + str(total_time) + '\n')
        log.close()

    if not args['disable_info']:
        if (args['num_clones'] > 0 or args['chrom_rate_super_clone']) and args['chrom_level_event']:
            record_clone_events(tree, os.path.join(args['out_path'], 'clone_events.tsv'), args['chrom_rate_super_clone'], clone_founders)
        record_events(tree, args['region_length'], os.path.join(args['out_path'], 'focal_events.tsv'))

if __name__ == '__main__':
    main()