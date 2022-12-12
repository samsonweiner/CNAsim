import os
import sys

class CustomError(Exception):
    pass

#Ensure inputs are sound
'''
def check_input(inputs):
    if out_path[-1] != '/':
        out_path += '/'
    if not os.path.isdir(out_path):
        call = subprocess.call(['mkdir', out_path])

    try:
        if tree_type not in [0, 1, 2]:
            raise CustomError
    except CustomError:
        print('Input Error: For tree method choose 0, 1, or 2')

    if tree_type == 2:
        try:
            if not os.path.exists(tree_path):
                raise CustomError
        except CustomError:
            print('Input Error: Tree file not found.')
'''

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
    return args

#gets chromosome and arm lengths from cytoband file
def hg38_chrom_lengths_from_cytoband(file_path, include_allosomes=False, include_arms=False):
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

# Summary sim stats 
def summary(tree, out_path, num_chroms, WGD, min_cn_length):
    f = open(out_path, 'w+')
    f.write('Summary\n-----------\n')
    f.write('WGD present: ' + str(WGD) + '\n')

    mutations = tree.get_mutations()

    whole_event_stats = {}
    CNV_count, CNV_lens, CNV_gains, CNV_copy = 0, 0, 0, 0
    for mut in mutations:
        if mut.category == 0:
            CNV_count += 1
            CNV_lens += mut.length
            CNV_gains += mut.event
            if mut.event == 1:
                CNV_copy += mut.copies
        elif mut.category == 1:
            if mut.cell not in whole_event_stats:
                whole_event_stats[mut.cell] = {}
            if mut.chrom not in whole_event_stats[mut.cell]:
                whole_event_stats[mut.cell][mut.chrom] = {}
            if mut.event == 0:
                whole_event_stats[mut.cell][mut.chrom][mut.allele] = [str(mut.cell), 'deletion', 'None']
            else:
                whole_event_stats[mut.cell][mut.chrom][mut.allele] = [str(mut.cell), 'gain', str(mut.copies)]

    f.write('\nWhole chromosomal events\n')
    for cell in whole_event_stats:
        for chrom in whole_event_stats[cell]:
            if 0 in whole_event_stats[cell][chrom]:
                f.write('Chrom' + str(chrom+1) + '-p: '+ '/'.join(whole_event_stats[cell][chrom][0]) + ' --- cell/event/copies' + '\n')
            if 1 in whole_event_stats[cell][chrom]:
                f.write('Chrom' + str(chrom+1) + '-m: '+ '/'.join(whole_event_stats[cell][chrom][1]) + ' --- cell/event/copies' + '\n')
    f.write('\nFocal events\n')
    if CNV_count != 0:
        f.write('Number of events: ' + str(CNV_count) + '\n')
        f.write('Average length: ' + str(round(CNV_lens / CNV_count, 3)) + ' * ' + str(min_cn_length) + '\n')
        f.write('Percent gain events: ' + str(round(CNV_gains / CNV_count, 3)) + '\n')
    if CNV_gains != 0:
        f.write('Average number of copies gained: ' + str(round(CNV_copy / CNV_gains, 3)) + '\n')
    
    '''
    tallies = {}
    for cell in tree.iter_leaves():
        total_unchanged, total_changed, num_changes, num_segments = 0, 0, 0, 0
        for chrom in range(num_chroms):
            for allele in [0, 1]:
                num_segments += len(cell.alterations[chrom][allele])
                for segment in cell.alterations[chrom][allele]:
                    if segment == 0:
                        total_unchanged += 1
                    if segment > 0:
                        total_changed += 1
                        num_changes += segment
        tallies[cell.name] = [total_unchanged / num_segments, total_changed / num_segments, num_changes / total_changed]

    global_averages = [0, 0, 0]
    for cell, stats in tallies.items():
        global_averages[0] += stats[0] 
        global_averages[1] += stats[1]
        global_averages[2] += stats[2]
    global_averages = [str(round(i / len(tallies), 3)) for i in global_averages]
    
    f.write('\nAlteration stats -- cell/fraction unchanged/fraction changed/average times changed\n')
    f.write('Average: ' + '/'.join(global_averages) + '\n')
    '''

    '''
    f.write('\nLost Chromosomes -- chrom/number of descendant leaves\n')
    found = {}
    for cell in tree.iter_descendants():
        for chrom in range(num_chroms):
            for allele in [0, 1]:
                if len(cell.genome[chrom][allele]) == 0:
                    k = str(chrom) + '_' + str(allele)
                    if k in found:
                        if not cell.name in found[k]:
                            f.write('Chrom' + k + ': ' + str(len(cell.iter_leaves())))
                            found[k] = [node.name for node in cell.iter_descendants()]
                    else:
                        f.write('Chrom' + k + ': ' + str(len(cell.iter_leaves())))
                        found[k] = [node.name for node in cell.iter_descendants()]
    '''
    f.close()

def record_cell_types(tree, out_path):
    f = open(out_path, 'w+')
    for n in tree.iter_descendants():
        f.write(n.name + '\t' + n.cell_type + '\n')
    f.close()

def record_events(tree, out_path):
    f = open(out_path, 'w+')
    f.write('\t'.join(['node', 'category', 'chrom', 'allele', 'homolog', 'arm', 'start', 'length', 'event', 'copies']) + '\n')
    for n in tree.iter_descendants():
        for e in n.events:
            f.write('\t'.join([str(x) for x in [n.name, e.category, e.chrom, e.allele, e.homolog, e.arm, e.start, e.length, e.event, e.copies]]) + '\n')
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