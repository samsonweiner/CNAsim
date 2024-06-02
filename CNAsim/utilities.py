import os
import sys

class InputError(Exception):
    def __init__(self, text):
        self.message = str(text)
        super().__init__(self.message)

class ChromNameError(Exception):
    def __init__(self):
        self.message = "Chromosome names in reference do not match the alternate and/or precomputed hg38 names. Check chromosome names or consider disabling the --use-hg38-static parameter with this reference."
        super().__init__(self.message)

class ChromNumError(Exception):
    def __init__(self):
        self.message = "Cannot have more than 22 chroms with the --use-hg38-static parameter toggled."
        super().__init__(self.message)

class ChromSizeError(Exception):
    def __init__(self):
        self.message = "Detected chromosomes in reference smaller than window size. Ensure the reference contains the desired chromosome sequences and remove alternate sequences if necessary."
        super().__init__(self.message)

class ModeError(Exception):
    def __init__(self):
        self.message = "The simulation mode needs to be specified correctly. Select 0 for CNP, 1 for seq, or 2 for both."
        super().__init__(self.message)

def convert(val):
    if val == 'True':
        return True
    if val == 'False':
        return False
    constructors = [int, float]
    for c in constructors:
        try:
            return c(val)
        except ValueError:
            pass
    return val

def handle_args(arguments):
    if arguments.param_file:
        args = {}
        f = open(arguments.param_file)
        lines = f.readlines()
        f.close()

        for line in lines:
            if line[0] == '[':
                param, val = line[line.index('(')+1:line.index(')')], line[line.index('[')+1:line.index(']')]
                args[param] = convert(val)

    else:
        args = vars(arguments)

    if not args['normal_fraction'] + args['pseudonormal_fraction'] <= 1.0:
        raise InputError('Normal and pseudornormal fractions cannot sum to more than 1.')
    return args

#gets chromosome and arm lengths from cytoband file
def hg38_chrom_lengths_from_cytoband_file(file_path, include_allosomes=False, include_arms=False):
    chrom_lens, arm_ratios = {}, {}

    f = open(file_path)
    lines = f.readlines()
    f.close

    for line in lines:
        info = line.strip().split('\t')
        chrom, start, end, arm = info[0], int(info[1]), int(info[2]), info[3][0]
        if chrom not in chrom_lens:
            chrom_lens[chrom] = {'p': 0, 'q': 0}
        chrom_lens[chrom][arm] = end

    if not include_allosomes:
        del chrom_lens['chrX']
        del chrom_lens['chrY']

    for chrom in chrom_lens:
        arm_ratios[chrom] = chrom_lens[chrom]['p'] / chrom_lens[chrom]['q']
        chrom_lens[chrom] = chrom_lens[chrom]['q']

    if include_arms:
        return chrom_lens, arm_ratios
    else:
        return chrom_lens
    
def hg38_chrom_lengths_from_cytoband(include_allosomes=False, include_arms=False):
    chrom_lens = {'chr1': 249250621, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr2': 243199373, 'chr20': 63025520, 'chr21': 48129895, 'chr22': 51304566, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chrX': 155270560, 'chrY': 59373566}
    arm_ratios = {'chr1': 0.5015, 'chr10': 0.2966, 'chr11': 0.39776, 'chr12': 0.26746, 'chr13': 0.15542, 'chr14': 0.16395, 'chr15': 0.18531, 'chr16': 0.40507, 'chr17': 0.29558, 'chr18': 0.22029, 'chr19': 0.44817, 'chr2': 0.38364, 'chr20': 0.43633, 'chr21': 0.27426, 'chr22': 0.28652, 'chr3': 0.45954, 'chr4': 0.26366, 'chr5': 0.26753, 'chr6': 0.35649, 'chr7': 0.3764, 'chr8': 0.31155, 'chr9': 0.34699, 'chrX': 0.39029, 'chrY': 0.21053}
    if not include_allosomes:
        del chrom_lens['chrX']
        del chrom_lens['chrY']
        del arm_ratios['chrX']
        del arm_ratios['chrY']
    if include_arms:
        return chrom_lens, arm_ratios
    else:
        return chrom_lens
    
def check_chrom_lengths(chrom_lens, window_size):
    for c,l in chrom_lens.items():
        if l < window_size:
            return False
    return True

# Summary sim stats 
def summary(tree, out_path, num_chroms, WGD, region_length):
    pass

def record_cell_types(tree, out_path, chrom_level_event, super_clone_rate, clone_founders=[]):
    f = open(out_path, 'w+')
    if super_clone_rate and chrom_level_event:
        #for c in tree.founder.children:
        #    if c not in clone_founders:
        #        clone_founders.append(c)
        clone_founders.extend(tree.founder.children)
    if len(clone_founders) > 0:
        clone_founders.sort(key=lambda x: len(x), reverse=True)
        cloneid = 0
        membership = {}
        for n in tree.iter_descendants():
            membership[n] = n.cell_type
        for clone in clone_founders:
            cloneid += 1
            for n in clone.iter_descendants():
                if n == clone:
                    membership[n] = 'clone' + str(cloneid) + '_founder'
                else:
                    membership[n] = 'clone' + str(cloneid)
        for n in tree.iter_descendants():
            f.write(n.name + '\t' + membership[n] + '\n')
    else:
        for n in tree.iter_descendants():
            f.write(n.name + '\t' + n.cell_type + '\n')
    f.close()
    return clone_founders

def record_clone_events(tree, out_path, super_clone_rate, clone_founders):
    f = open(out_path, 'w+')
    cloneids = {clone_founders[i]: 'clone' + str(i+1) for i in range(len(clone_founders))}
    eventdict = {0: 'del', 1: 'dup'}
    scaledict = {1: {None: 'whole'}, 2: {'p': 'p', 'q': 'q', None: 'Other'}}
    f.write('\t'.join(['cellname', 'cloneid', 'chrom', 'allele', 'scale', 'event']) + '\n')
    for e in tree.founder.events:
        if e.category != 0:
            f.write('\t'.join([str(x) for x in ['founder', 'None', e.chrom, e.allele, scaledict[e.category][e.arm], eventdict[e.event]]]) + '\n')
    for c in clone_founders:
        for e in c.events:
            if e.category != 0:
                f.write('\t'.join([str(x) for x in [c.name, cloneids[c], e.chrom, e.allele, scaledict[e.category][e.arm], eventdict[e.event]]]) + '\n')
    f.close()

def record_events(tree, region_length, out_path):
    f = open(out_path, 'w+')
    eventdict = {0: 'del', 1: 'gain'}
    f.write('\t'.join(['cellname', 'chrom', 'allele', 'tumor_start_idx', 'ref_start_idx', 'length', 'event', 'copies']) + '\n')
    for n in tree.iter_descendants():
        for e in n.events:
            if e.category == 0:
                f.write('\t'.join([str(x) for x in [n.name, e.chrom, e.allele, e.start*region_length + 1, e.ref_idx*region_length + 1, e.length*region_length, eventdict[e.event], e.copies]]) + '\n')
    f.close()

def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size