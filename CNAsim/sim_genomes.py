import numpy as np
from .tree import Tree
from collections import Counter
import os
import pickle


# Categories --> 0: focal, 1: whole-chromosomal, 2: chromosome-arm
class CNV():
    def __init__(self, cell = None, category=None, chrom=None, allele=None, homolog=None, arm=None, start=None, length=None, event=None, copies=None, ref_idx=None):
        self.cell = cell
        self.category = category
        self.chrom = chrom
        self.allele = allele
        self.start = start
        self.length = length
        self.event = event
        self.copies = copies
        self.homolog = homolog
        self.arm = arm
        self.ref_idx = ref_idx

def get_chrom_proportions(sequence):
    combined_chrom_len = dict(zip(sequence.keys(), map(lambda x: len(sequence[x][0]) + len(sequence[x][1]), sequence.keys())))
    total_len = sum(combined_chrom_len.values())
    combined_chrom_proportions = dict(zip(combined_chrom_len.keys(), [val/total_len for val in combined_chrom_len.values()]))
    return combined_chrom_proportions

# Creates initial genome from either from reference file, hg38 lengths, or with fixed chrom lengths and arm ratios
def init_diploid_genome(region_length, chrom_names, chrom_lens, arm_ratios):
    genome, regions = {}, {}
    static1, static2 = isinstance(chrom_lens, dict), isinstance(arm_ratios, dict)
    for chrom in chrom_names:
        if static1:
            chrom_len = chrom_lens[chrom]
        else:
            chrom_len = chrom_lens
        if static2:
            arm_ratio = arm_ratios[chrom]
        else:
            arm_ratio = arm_ratios

        num_regions = round(chrom_len / region_length)
        cent_idx = round(num_regions*arm_ratio)
        hap_profile = [i for i in range(num_regions)]
        hap_profile.insert(cent_idx, 'X')
        chrom_profile = [[[i for i in hap_profile]], [[i for i in hap_profile]]]
        genome[chrom] = chrom_profile
        regions[chrom] = num_regions
    return genome, regions



# Alters the genome of a cell based on the assignmened whole chrom events.
#def gen_whole_chrom_events(cell, chrom_names, chrom_arm_loc, chrom_event_rate, chrom_arm_rate, chrom_event_type, chrom_event_copy):
#    for chrom in chrom_names:
#        for allele in [0, 1]:
#            if np.random.binomial(1, chrom_event_rate) == 1:
#                if np.random.binomial(1, chrom_arm_rate) == 1:
#                    num_regions = len(cell.genome[chrom][allele])
#                    arm_region = chrom_arm_loc[chrom]
#                    while arm_region not in cell.genome[chrom][allele]:
#                        pass
#                    arm_idx = num_regions - cell.genome[chrom][allele][::-1].index(chrom_arm_loc[chrom])
#                    event_type = np.random.binomial(1, chrom_event_type)
#                    if np.random.binomial(1, 0.5) == 1:
#                        if event_type == 0:
#                            cell.genome[chrom][allele] = cell.genome[chrom][allele][:arm_idx]
#                    else:
#                        start_idx = arm_idx
#                        end_idx = num_regions
#                    event_type = np.random.binomial(1, chrom_event_type)
#                else:
#                    event_type = np.random.binomial(1, chrom_event_type)
#                    if event_type == 0:
#                        cell.genome[chrom][allele] = []
#                        cell.events.append(CNV(cell=cell, category=1, chrom=chrom, allele=allele, event=event_type))
#                    else:
#                        num_copies = np.random.geometric(chrom_event_copy)
#                        cell.genome[chrom][allele] = num_copies*cell.genome[chrom][allele]
#                        cell.events.append(CNV(cell=cell, category=1, chrom=chrom, allele=allele, event=event_type, copies=num_copies))

def gen_whole_chrom_event(cell, chrom_names, chrom_event_type):
    intact_chroms = []
    for chrom in chrom_names:
        for allele in [0, 1]:
            # Only consider chromosomes with both arms intact. If no X is present, can assume one of the arms was lost already.
            intact_chroms.extend([(chrom, allele, h) for h in range(len(cell.genome[chrom][allele])) if 'X' in cell.genome[chrom][allele][h]])
    if len(intact_chroms) == 0:
        print('No intact chromosomes remaining.')
        return
    
    chrom_idx = np.random.randint(len(intact_chroms))
    chrom, allele, homolog = intact_chroms[chrom_idx]

    event_type = np.random.binomial(1, chrom_event_type)
    if event_type == 0:
        del cell.genome[chrom][allele][homolog]
    else:
        cell.genome[chrom][allele].append([i for i in cell.genome[chrom][allele][homolog]])
    
    cell.events.append(CNV(cell=cell, category=1, chrom=chrom, allele=allele, homolog=homolog, event=event_type))

def gen_chrom_arm_event(cell, chrom_names, chrom_event_type):
    intact_arms = []
    for chrom in chrom_names:
        for allele in [0, 1]:
            for h in range(len(cell.genome[chrom][allele])):
                if len(cell.genome[chrom][allele][h]) > 0:
                    if 'X' in cell.genome[chrom][allele][h]:
                        split_idx = cell.genome[chrom][allele][h].index('X')
                        #if split_idx != 0:
                        intact_arms.append((chrom, allele, h, 'p'))
                        #if split_idx != len(cell.genome[chrom][allele][h]) - 1:
                        intact_arms.append((chrom, allele, h, 'q'))
                    else:
                        # Arm must have already been removed, but still consider the other one
                        intact_arms.append((chrom, allele, h, None))
    if len(intact_arms) == 0:
        print('No intact arms remaining.')
        return

    arm_idx = np.random.randint(len(intact_arms))
    chrom, allele, homolog, arm = intact_arms[arm_idx]
    if arm != None:
        split_idx = cell.genome[chrom][allele][homolog].index('X')
    
    event_type = np.random.binomial(1, chrom_event_type)
    if event_type == 0:
        if arm == 'p':
            cell.genome[chrom][allele][homolog] = cell.genome[chrom][allele][homolog][split_idx+1:]
        elif arm == 'q':
            cell.genome[chrom][allele][homolog] = cell.genome[chrom][allele][homolog][:split_idx]
        else:
            del cell.genome[chrom][allele][homolog]
    else:
        if arm == 'p':
            cell.genome[chrom][allele][homolog] = 2 * cell.genome[chrom][allele][homolog][:split_idx] + cell.genome[chrom][allele][homolog][split_idx:]
        elif arm == 'q':
            cell.genome[chrom][allele][homolog] = cell.genome[chrom][allele][homolog][:split_idx + 1] + 2 * cell.genome[chrom][allele][homolog][split_idx+1:]
        else:
            cell.genome[chrom][allele][homolog] = 2 * cell.genome[chrom][allele][homolog]
    
    cell.events.append(CNV(cell=cell, category=2, chrom=chrom, allele=allele, homolog=homolog, arm=arm, event=event_type))

def gen_focal_event(cell, chrom_names, length_mean, min_length, event_rate, copy_param):
    #Determine chrom and allele of event from those that haven't been lost.
    #CN_allele = np.random.binomial(1, 0.5)
    #chrom_lens = [len(cell.genome[chrom][CN_allele]) for chrom in chrom_names]
    intact_chroms = []
    for chrom in chrom_names:
            for allele in [0, 1]:
                #print(cell.name, chrom, allele, cell.genome[chrom][allele])
                intact_chroms.extend([(chrom, allele, h) for h in range(len(cell.genome[chrom][allele])) if len(cell.genome[chrom][allele][h]) > 0])
    chrom_proportions = [len(cell.genome[x[0]][x[1]][x[2]]) - 1 for x in intact_chroms]
    chrom_proportions = [x/sum(chrom_proportions) for x in chrom_proportions]
    chrom_idx = np.random.choice(list(range(len(intact_chroms))), p=chrom_proportions)
    CN_chrom, CN_allele, CN_homolog = intact_chroms[chrom_idx]

    #if sum(chrom_lens) == 0:
    #    print('Somehow all chroms of allele ' + str(CN_allele) + ' have been lost.')
    #chrom_proportions = [i / sum(chrom_lens) for i in chrom_lens]
    #CN_chrom = np.random.choice(chrom_names, p=chrom_proportions)

    #Determine properties of event
    num_regions = len(cell.genome[CN_chrom][CN_allele][CN_homolog])
    CN_size = max(min(round(np.random.exponential(length_mean)), num_regions), min_length)
    CN_type = np.random.binomial(1, event_rate)
    CN_copies = np.random.geometric(copy_param)

    #get starting location
    start_idx = np.random.randint(num_regions - CN_size + 1)
    ref_idx = cell.genome[CN_chrom][CN_allele][CN_homolog][start_idx]
    end_idx = start_idx + CN_size
    
    update_arm = False
    if 'X' in cell.genome[CN_chrom][CN_allele][CN_homolog]:
        arm_idx = cell.genome[CN_chrom][CN_allele][CN_homolog].index('X')
        if arm_idx >= start_idx and arm_idx <= end_idx:
            update_arm = True

    #update genome
    if CN_type == 1:
        gain_segment = CN_copies * cell.genome[CN_chrom][CN_allele][CN_homolog][start_idx:end_idx]
        if update_arm:
            gain_segment = list(filter(lambda x: x != 'X', gain_segment))
        cell.genome[CN_chrom][CN_allele][CN_homolog] = cell.genome[CN_chrom][CN_allele][CN_homolog][:end_idx] + gain_segment + cell.genome[CN_chrom][CN_allele][CN_homolog][end_idx:]
    else:
        cell.genome[CN_chrom][CN_allele][CN_homolog] = cell.genome[CN_chrom][CN_allele][CN_homolog][:start_idx] + cell.genome[CN_chrom][CN_allele][CN_homolog][end_idx:]
        if update_arm and CN_size != num_regions:
            cell.genome[CN_chrom][CN_allele][CN_homolog].insert(start_idx, 'X')
        if len(cell.genome[CN_chrom][CN_allele][CN_homolog]) == 0:
            del cell.genome[CN_chrom][CN_allele][CN_homolog]
        else:
            if cell.genome[CN_chrom][CN_allele][CN_homolog][0] == 'X':
                cell.genome[CN_chrom][CN_allele][CN_homolog].pop(0)
            if cell.genome[CN_chrom][CN_allele][CN_homolog][-1] == 'X':
                cell.genome[CN_chrom][CN_allele][CN_homolog].pop(-1)

    cell.events.append(CNV(cell=cell, category=0, chrom=CN_chrom, allele=CN_allele, start=start_idx, length=CN_size, event=CN_type, copies=CN_copies, ref_idx=ref_idx))

def mutate_genome(node, args, chrom_names):
    # Add WGD to founder node, if necessary
    if node.cell_type == 'founder':
        num_events = len(node.events)
        node.events = []
        if args['WGD']:
            for chrom in chrom_names:
                for allele in [0, 1]:
                    node.genome[chrom][allele].append([x for x in node.genome[chrom][allele][0]])

    # Add whole chrom events to founder and it's children, if necessary.
    if node.cell_type == 'founder' or node.parent.cell_type == 'founder' or node.cell_type == 'clone':
        if args['chrom_level_event']:
            if node.cell_type == 'clone':
                num_chrom_events = max(np.random.poisson(args['chrom_rate_clone']), 1) # we want to make sure theres atleast 1 event
            elif node.cell_type == 'founder':
                num_chrom_events = np.random.poisson(args['chrom_rate_founder']) # --> Changing implementation where this is a param specific to founder edges
            else:
                num_chrom_events = np.random.poisson(args['chrom_rate_super_clone'])
            for i in range(num_chrom_events):
                chrom_alt_type = np.random.binomial(1, args['chrom_arm_rate'])
                if chrom_alt_type == 1: #arm-level event
                    gen_chrom_arm_event(node, chrom_names, args['chrom_event_type'])
                else: #whome-chrom event
                    gen_whole_chrom_event(node, chrom_names, args['chrom_event_type'])
            #gen_whole_chrom_events(node, chrom_names, args['whole_chrom_rate'], args['whole_chrom_type'], args['whole_chrom_copy'])
    
    # Generates focal events. If above aneuploid tree, number of events precomputed.
    if node.cell_type == 'pseudonormal':
        num_events = len(node.events)
        node.events = []
    elif node.cell_type != 'founder': #foudner events already determined
        if args['placement_type'] == 0:
            num_events = np.random.poisson(args['placement_param'])
        if args['placement_type'] == 1:
            num_events = np.random.poisson(node.length)
        if args['placement_type'] == 2:
            num_events = int(args['placement_param'])
    #if node.is_root():
        #num_events = int(num_events * args['root_event_mult'])
    for i in range(num_events):
        gen_focal_event(node, chrom_names, args['cn_length_mean']/args['region_length'], max(round(args['min_cn_length']/args['region_length']), 1), args['cn_event_rate'], args['cn_copy_param'])

def genome_to_profile(node, chrom_names, num_regions, bins):
    # Collapsing genome in the form of a counter
    for chrom in chrom_names:
        for allele in [0, 1]:
            node.genome[chrom][allele] = [x for homolog in node.genome[chrom][allele] for x in homolog if x != 'X']
            node.genome[chrom][allele] = dict(Counter(node.genome[chrom][allele]))
            for i in set([i for i in range(num_regions[chrom])]) - set(node.genome[chrom][allele].keys()):
                node.genome[chrom][allele][i] = 0
    
    # Find the average number of copies over the regions in each bin
    node.profile = {}
    for chrom in chrom_names:
        node.profile[chrom] = [[], []]
        for allele in [0, 1]:
            for i in range(len(bins[chrom])-1):
                regions = list(range(bins[chrom][i], bins[chrom][i+1]))
                node.profile[chrom][allele].append(round(sum([node.genome[chrom][allele][j] for j in regions]) / len(regions)))

def evolve_tree(node, args, chrom_names, num_regions, bins=None):
    if not node.is_root():
        node.inheret()
    
    if not node.cell_type == 'normal':
        mutate_genome(node, args, chrom_names)
    
    if args['mode'] == 1 or args['mode'] == 2:
        with open(os.path.join(args['out_path'], node.name + '.pkl'), 'wb') as f:
            pickle.dump(node.genome, f)

    #if (args['mode'] == 0 or args['mode'] == 2) and bins != None:
    if node.is_leaf():
        genome_to_profile(node, chrom_names, num_regions, bins)
    else:
        with open(os.path.join(args['out_path'], node.name + '.pkl'), 'wb') as f:
            pickle.dump(node.genome, f)
            
    for child in node.children:
        evolve_tree(child, args, chrom_names, num_regions, bins=bins)

    del node.genome

def prepare_ancestral_profiles(tree, args, chrom_names, num_regions, bins=None):
    for node in tree.iter_preorder():
        if not node.is_leaf():
            with open(os.path.join(args['out_path'], node.name + '.pkl'), 'rb') as f:
                node.genome = pickle.load(f)
            genome_to_profile(node, chrom_names, num_regions, bins)
            del node.genome
            os.remove(os.path.join(args['out_path'], node.name + '.pkl'))



