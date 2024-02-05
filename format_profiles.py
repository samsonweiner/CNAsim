import numpy as np
from collections import Counter

def get_genomes(tree):
    genomes = {}
    for leaf in tree.iter_leaves():
        genomes[leaf.name] = leaf.genome
    return genomes

def collapse_genomes(genomes, normal_diploid_genome, num_chroms):
    for chrom in range(num_chroms):
        for leaf, genome in genomes.items():
            for allele in [0, 1]:
                genome[chrom][allele] = dict(Counter(genome[chrom][allele]))
                for i in set(normal_diploid_genome[chrom][allele]) - set(genome[chrom][allele].keys()):
                    genome[chrom][allele][i] = 0
    return genomes
                
def format_CN_profiles(tree, normal_diploid_genome, num_chroms, region_length, bin_len):
    # Consolidate leaf genomes
    genomes = get_genomes(tree)
    genome_counters = collapse_genomes(genomes, normal_diploid_genome, num_chroms)

    # Assign regions to bins
    bins, bin_coords = {}, {}
    regions_per_bin = np.floor(bin_len / region_length)
    for chrom in range(num_chroms):
        bins[chrom], bin_coords[chrom] = {}, {}
        bin_count = 0
        for i in normal_diploid_genome[chrom][0]:
            if i == 0:
                bins[chrom][bin_count] = [i]
            else:
                if i % regions_per_bin == 0:
                    bin_coords[chrom][bin_count] = [bins[chrom][bin_count][0]*region_length, (bins[chrom][bin_count][-1]+1)*region_length]
                    bin_count += 1
                    bins[chrom][bin_count] = [i]
                else:
                    bins[chrom][bin_count].append(i)
        bin_coords[chrom][bin_count] = [bins[chrom][bin_count][0]*region_length, (bins[chrom][bin_count][-1]+1)*region_length]

    # Find the average number of copies over the regions in each bin
    CN_profiles = {}
    for chrom in range(num_chroms):
        CN_profiles[chrom] = {}
        sorted_bins = list(bins[chrom].keys())
        sorted_bins.sort()
        for leaf, genome in genome_counters.items():
            CN_profiles[chrom][leaf] = [[], []]
            for b in sorted_bins:
                regions = bins[chrom][b]
                for allele in [0, 1]:
                    CN_profiles[chrom][leaf][allele].append(round(sum([genome[chrom][allele][i] for i in regions]) / len(regions)))
    
    return CN_profiles, bin_coords

def save_CN_profiles_leaves(tree, chrom_names, bins, region_length, filepath):
    f = open(filepath, 'w+')
    headers = ['CELL', 'chrom', 'start', 'end', 'CN states']
    f.write('\t'.join(headers) + '\n')

    leaves = [leaf for leaf in tree.iter_leaves()]

    # Writes to file in order of chroms, bins, cell
    for chrom in chrom_names:
        num_bins = len(bins[chrom]) - 1
        for i in range(num_bins):
            bin_start, bin_end = str(bins[chrom][i]*region_length), str(bins[chrom][i+1]*region_length)
            for leaf in leaves:
                CN_state = str(leaf.profile[chrom][0][i]) + ',' + str(leaf.profile[chrom][1][i])
                line = [leaf.name, chrom, bin_start, bin_end, CN_state]
                f.write('\t'.join(line) + '\n')
    f.close()

def save_CN_profiles_ancestors(tree, chrom_names, bins, region_length, filepath):
    f = open(filepath, 'w+')
    headers = ['CELL', 'chrom', 'start', 'end', 'CN states']
    f.write('\t'.join(headers) + '\n')

    ancestors = [node for node in tree.iter_preorder() if not node.is_leaf()]

    # Writes to file in order of chroms, bins, cell
    for chrom in chrom_names:
        num_bins = len(bins[chrom]) - 1
        for i in range(num_bins):
            bin_start, bin_end = str(bins[chrom][i]*region_length), str(bins[chrom][i+1]*region_length)
            for ancestor in ancestors:
                CN_state = str(ancestor.profile[chrom][0][i]) + ',' + str(ancestor.profile[chrom][1][i])
                line = [ancestor.name, chrom, bin_start, bin_end, CN_state]
                f.write('\t'.join(line) + '\n')
    f.close()
            
