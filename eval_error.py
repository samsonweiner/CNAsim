import argparse
import os

def process_data(data_path):
    f = open(data_path)
    lines = f.readlines()
    f.close()

    data, cell_names, chrom_names, bin_names = {}, [], [], {}
    profiles = {}
    for row in lines[1:]:
        cell, chrom, start, end, CN = row[:-1].split('\t')
        start = int(start)
        cn_a, cn_b = int(CN[1:CN.index(',')]), int(CN[CN.index(',')+1:-1])
        if cell not in data:
            data[cell] = {}
            cell_names.append(cell)
        if chrom not in data[cell]:
            data[cell][chrom] = {0: {}, 1: {}}
            profiles[chrom] = {}
            chrom_names.append(chrom)
        if chrom not in bin_names:
            bin_names[chrom] = set()
        bin_names[chrom].add(start)
        data[cell][chrom][0][start], data[cell][chrom][1][start] = cn_a, cn_b

    for chrom in chrom_names:
        bin_names[chrom] = list(bin_names[chrom])
        bin_names[chrom].sort()
        for allele in [0, 1]:
            profiles[chrom][allele] = [[data[cell][chrom][allele][b] for b in bin_names[chrom]] for cell in cell_names]
            
    return profiles, cell_names, bin_names


def profile_noise_eval(noisy_prof, clean_prof, cell_names, tolerance):
    gt_bkpt, ns_bkpt = {}, {}
    for chrom in clean_prof:
        gt_bkpt[chrom] = {0: {}, 1: {}}
        ns_bkpt[chrom] = {0: {}, 1: {}}
        for allele in [0, 1]:
            for cell in range(len(cell_names)):
                cur_gt_profile = clean_prof[chrom][allele][cell]
                cur_ns_profile = noisy_prof[chrom][allele][cell]
                gt_bkpt[chrom][allele][cell] = [0 for b in range(len(cur_gt_profile) - 1)]
                ns_bkpt[chrom][allele][cell] = [0 for b in range(len(cur_gt_profile) - 1)]
                for b in range(len(cur_gt_profile) - 1):
                    if cur_gt_profile[b] != cur_gt_profile[b+1]:
                        gt_bkpt[chrom][allele][cell][b] = cur_gt_profile[b+1] - cur_gt_profile[b]
                    if cur_ns_profile[b] != cur_ns_profile[b+1]:
                        ns_bkpt[chrom][allele][cell][b] = cur_ns_profile[b+1] - cur_ns_profile[b]
                

    TP, FP, FN = 0, 0, 0
    for chrom in ns_bkpt:
        for allele in [0, 1]:
            for cell in range(len(cell_names)):
                for loc, change in enumerate(ns_bkpt[chrom][allele][cell]):
                    if change != 0:
                        pool = [loc]
                        if loc > (-1 + tolerance) and tolerance > 0:
                            pool.append(loc - 1)
                        if loc < len(ns_bkpt[chrom][allele][cell]) - (tolerance) and tolerance > 0:
                            pool.append(loc + 1)
                        
                        match = False
                        for candidate in pool:
                            if (change > 0 and gt_bkpt[chrom][allele][cell][candidate] > 0) or (change < 0 and gt_bkpt[chrom][allele][cell][candidate] < 0):
                                match = True
                                gt_bkpt[chrom][allele][cell][candidate] = 0
                                break
                        if match:
                            TP += 1
                        else:
                            FP += 1
    
    for chrom in gt_bkpt:
        for allele in [0, 1]:
            for cell in range(len(cell_names)):
                FN += sum([1 for i in gt_bkpt[chrom][allele][cell] if i != 0])
    
    recall = round(TP / (TP + FN), 3)
    precision = round(TP / (TP + FP), 3)

    return recall, precision


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', type=str, required=True, help='Path to dataset directory.')
parser.add_argument('-t', '--tolerance', type=int, default=0, help='Tolerance for accepting a candidate breakpoint as correct in terms of number of bins from correct location.')
args = parser.parse_args()

noisy_prof, c1, b1 = process_data(os.path.join(args.path, 'profiles.tsv'))
clean_prof, c2, b2 = process_data(os.path.join(args.path, 'clean_profiles.tsv'))
recall, precision = profile_noise_eval(noisy_prof, clean_prof, c1, args.tolerance)

print('Recall:', recall)
print('Precision:', precision)