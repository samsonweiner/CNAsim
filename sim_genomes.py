import numpy as np
from tree import Tree
import copy
from collections import Counter

from utilities import hg38_chrom_lengths_from_cytoband

# Categories --> 0: focal, 1: chromosomal, 2: WGD
class CNV():
    def __init__(self, cell = None, category=None, chrom=None, allele=None, start=None, length=None, event=None, copies=None):
        self.cell = cell
        self.category = category
        self.chrom = chrom
        self.allele = allele
        self.start = start
        self.length = length
        self.event = event
        self.copies = copies

def get_chrom_proportions(sequence):
    combined_chrom_len = dict(zip(sequence.keys(), map(lambda x: len(sequence[x][0]) + len(sequence[x][1]), sequence.keys())))
    total_len = sum(combined_chrom_len.values())
    combined_chrom_proportions = dict(zip(combined_chrom_len.keys(), [val/total_len for val in combined_chrom_len.values()]))
    return combined_chrom_proportions

def scale_edge_lengths(tree, place_param):
    leaf_edge_lens = [leaf.length for leaf in tree.iter_leaves()]
    avg_leaf_len = sum(leaf_edge_lens) / len(leaf_edge_lens)
    scalar = place_param / avg_leaf_len

    for node in tree.iter_descendants():
        if node.is_root():
            node.length = place_param
        else:
            node.length = node.length * scalar

def init_diploid_genome(min_cn_len, chrom_names, chrom_len, use_hg38):
    genome = {}
    if use_hg38:
        chrom_lens = hg38_chrom_lengths_from_cytoband('resources/cytoBand.txt', include_allosomes=False, include_arms=False)
    for chrom in chrom_names:
        if use_hg38:
            num_regions = round(chrom_lens[chrom] / min_cn_len)
        else:
            num_regions = round(chrom_len / min_cn_len)
        chrom_profile = [[i for i in range(num_regions)], [i for i in range(num_regions)]]
        genome[chrom] = chrom_profile
    return genome

# Alters the genome of a cell based on the assignmened whole chrom events.
def gen_whole_chrom_events(cell, chrom_names, whole_chrom_rate, whole_chrom_type, whole_chrom_copy):
    for chrom in chrom_names:
        for allele in [0, 1]:
            if np.random.binomial(1, whole_chrom_rate) == 1:
                event_type = np.random.binomial(1, whole_chrom_type)
                num_copies = np.random.geometric(whole_chrom_copy)

                if event_type == 0:
                    cell.genome[chrom][allele] = []
                    cell.events.append(CNV(cell=cell, category=1, chrom=chrom, allele=allele, event=event_type))
                else:
                    cell.genome[chrom][allele] = num_copies*cell.genome[chrom][allele]
                    cell.events.append(CNV(cell=cell, category=1, chrom=chrom, allele=allele, event=event_type, copies=num_copies))


def gen_focal_event(cell, chrom_names, length_mean, event_rate, copy_param):
    #Determine chrom and allele of event from those that haven't been lost.
    CN_allele = np.random.binomial(1, 0.5)
    chrom_lens = [len(cell.genome[chrom][CN_allele]) for chrom in chrom_names]
    if sum(chrom_lens) == 0:
        print('Somehow all chroms of allele ' + str(CN_allele) + ' have been lost.')
    else:
        chrom_proportions = [i / sum(chrom_lens) for i in chrom_lens]
        CN_chrom = np.random.choice(chrom_names, p=chrom_proportions)

        #Determine properties of event
        num_regions = len(cell.genome[CN_chrom][CN_allele])
        CN_size = min(round(np.random.exponential(length_mean)), round(0.5*num_regions))
        CN_type = np.random.binomial(1, event_rate)
        CN_copies = np.random.geometric(copy_param)

        #get starting location
        start_idx = np.random.randint(num_regions - CN_size + 1)
        end_idx = start_idx + CN_size

        #update genome
        if CN_type == 1:
            gain_segment = CN_copies * cell.genome[CN_chrom][CN_allele][start_idx:end_idx]
            cell.genome[CN_chrom][CN_allele] = cell.genome[CN_chrom][CN_allele][:end_idx] + gain_segment + cell.genome[CN_chrom][CN_allele][end_idx:]
        else:
            cell.genome[CN_chrom][CN_allele] = cell.genome[CN_chrom][CN_allele][:start_idx] + cell.genome[CN_chrom][CN_allele][end_idx:]

        cell.events.append(CNV(cell=cell, category=0, chrom=CN_chrom, allele=CN_allele, start=start_idx, length=CN_size, event=CN_type, copies=CN_copies))

def mutate_genome(node, args, chrom_names):
    # Add whole chrom events to root and it's children, if necessary.
    if node.is_root() or node.parent.is_root():
        if args['whole_chrom_event']:
                gen_whole_chrom_events(node, chrom_names, args['whole_chrom_rate'], args['whole_chrom_type'], args['whole_chrom_copy'])
    
    # Generates focal events. If method == 1, proportions are precomputed
    if args['placement_type'] == 0:
        num_events = np.random.poisson(args['placement_param'])
    if args['placement_type'] == 1:
        num_events = np.random.poisson(node.length)
    if args['placement_type'] == 2:
        num_events = int(args['placement_param'])
        
    if node.is_root():
        num_events = int(num_events * args['root_event_mult'])

    for i in range(num_events):
        gen_focal_event(node, chrom_names, args['cn_length_mean']/args['min_cn_length'], args['cn_event_rate'], args['cn_copy_param'])


def format_profile(node, chrom_names, init_genome, bins):
    # Collapsing genome in the form of a counter
    for chrom in chrom_names:
        for allele in [0, 1]:
            node.genome[chrom][allele] = dict(Counter(node.genome[chrom][allele]))
            for i in set(init_genome[chrom][allele]) - set(node.genome[chrom][allele].keys()):
                node.genome[chrom][allele][i] = 0
    
    # Find the average number of copies over the regions in each bin
    node.profile = {}
    for chrom in chrom_names:
        node.profile[chrom] = [[], []]
        for allele in [0, 1]:
            for i in range(len(bins[chrom])-1):
                regions = list(range(bins[chrom][i], bins[chrom][i+1]))
                node.profile[chrom][allele].append(round(sum([node.genome[chrom][allele][i] for i in regions]) / len(regions)))


def evolve_tree(node, args, chrom_names, init_genome, bins):
    if not node.is_root():
        node.inheret()

    mutate_genome(node, args, chrom_names)

    if node.is_leaf():
        format_profile(node, chrom_names, init_genome, bins)
    else:
        evolve_tree(node.children[0], args, chrom_names, init_genome, bins)
        evolve_tree(node.children[1], args, chrom_names, init_genome, bins)

    del node.genome

#master function
def gen_profiles(args, tree, chrom_names, init_genome):
    # Initialization
    if args['placement_type'] == 1:
        scale_edge_lengths(tree, args['placement_param'])

    tree.root.genome = copy.deepcopy(init_genome)
    chrom_names = ['chr' + str(i+1) for i in range(args['num_chroms'])]
    if args['WGD']:
        for chrom in chrom_names:
            for allele in [0, 1]:
                tree.root.genome[chrom][allele] = 2 * tree.root.genome[chrom][allele]
    
    # Assign regions to bins
    regions_per_bin = np.floor(args['bin_length']/args['min_cn_length'])
    bins = {}
    for chrom in chrom_names:
        bins[chrom] = [k for k in init_genome[chrom][0] if k % regions_per_bin == 0]
        bins[chrom][-1] = init_genome[chrom][0][-1]

    # begin evolution from root
    evolve_tree(tree.root, args, chrom_names, init_genome, bins)

    return bins

