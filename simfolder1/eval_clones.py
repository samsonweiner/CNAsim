from willowtree import *
import argparse
import os

output_paths = {
    'medicc2': ['medicc2/medicc_input_final_tree.new', 'medicc2/times.txt', 0],
    'medalt': ['MEDALT/lineage_tree.nwk', 'MEDALT/times.txt', 0],
    'cnp2cnp': ['cnp2cnp/tree.nwk', 'cnp2cnp/times.txt', 0],
    'standard-euc': ['standard/euclidean_tree.nwk', 'standard/times.txt', 0],
    'standard-man': ['standard/manhattan_tree.nwk', 'standard/times.txt', 1],
    'breakpoint-man': ['breakpoint/manhattan_tree.nwk', 'breakpoint/times.txt', 0],
    'breakpoint-root': ['breakpoint/root_tree.nwk', 'breakpoint/times.txt', 1],
    'breakpoint-log': ['breakpoint/log_tree.nwk', 'breakpoint/times.txt', 2],
}

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--option', action='store_true')
parser.add_argument('-i', '--input', type=str, required=True, help='Path to dataset.')
args = parser.parse_args()

def validate(data_path):
    log = open(data_path + 'clone_results.txt', 'w+')

    gt = Tree(newick=data_path + 'tree.nwk')
    node_types = open(data_path+'cell_types.tsv')
    lines = node_types.readlines()
    node_types.close()
    clone_names = []
    for line in lines:
        info = line.strip().split('\t')
        if info[1] == 'clone':
            clone_names.append(info[0])
    
    clones = []
    for c in clone_names:
        n = gt.ref_node(c)
        clones.append([leaf.name for leaf in n.get_leaves()])
    
    for method, paths in output_paths.items():
        if os.path.exists(data_path + paths[0]):
            rec_tree = Tree(newick=data_path + paths[0])
            score = avg_f1_tree(rec_tree, clones)
            log.write('\t'.join([method, str(score)]) + '\n')

def compile(data_path):
    results = {}
    for i in range(1, 21):
        f = open(data_path + str(i) + '/clone_results.txt')
        lines = f.readlines()
        f.close()

        for line in lines:
            info = line.split('\t')
            method = info[0]
            if method not in results:
                results[method] = []
            results[method].append(float(info[1]))
    
    print(data_path)
    for method, res in results.items():
        avg_score = round(sum(results[method]) / len(results[method]), 3)
        print(method + ': ', avg_score)

def compute_f1(node, clone):
    leaves = [leaf.name for leaf in node.get_leaves()]
    TP, FP, FN = 0, 0, 0
    for leaf in leaves:
        if leaf in clone:
            TP += 1
        else:
            FP += 1
    for c in clone:
        if c not in leaves:
            FN += 1

    f1 = TP / (TP + 0.5*(FP + FN))
    return f1

def find_best_match(tree, clone):
    cur_best = None
    cur_score = 0

    for n in tree.iter_descendants():
        if not n.is_leaf() and not n.is_root():
            s = compute_f1(n, clone)
            if s > cur_score:
                cur_best = n
                cur_score = s
    
    return cur_best, cur_score

def avg_f1_tree(tree, clones):
    scores = []
    for clone in clones:
        n, s = find_best_match(tree, clone)
        scores.append(s)
    return sum(scores) / len(scores)

if args.input[-1] != '/':
    args.input += '/'

if args.option:
    compile(args.input)
else:
    validate(args.input)