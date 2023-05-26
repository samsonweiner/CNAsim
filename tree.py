import numpy as np
import scipy.stats
import subprocess
import msprime as ms
import os
from collections import deque, Counter
import copy

class Node:
    def __init__(self, name='', edge_len = 0, parent = None, cell_type='aneuploid'):
        self.name = name
        self.length = edge_len
        self.parent = parent
        self.children = []
        self.events = []
        self.cell_type = cell_type
        #self.whole_chrom_events = []
        self.genome = None
        #self.alterations = None
        self.profile = None

    def __str__(self):
        if self.name:
            return str(self.name)
        else:
            return ''

    def __len__(self):
        return len(self.get_leaves())
    
    #format codes --> 0: leaf names only, 1: leaf names + lengths only, 2: leaf and internal names, 3: leaf and internal names + lengths
    def write_newick(self, terminate=True, format=0):
        if self.is_leaf():
            if format == 0 or format == 2:
                return self.name
            elif format == 1 or format == 3:
                return self.name + ':' + str(self.length)
        else:
            newick_str = '('
            for child in self.children:
                newick_str += child.write_newick(terminate=False, format=format) + ','
            newick_str = newick_str[:-1]
            newick_str += ')'
            if format == 2 or format == 3 or terminate:
                newick_str += self.name
            if (format == 1 or format == 3) and not terminate :
                newick_str += ':' + str(self.length)
        if terminate:
            return newick_str + ';'
        else:
            return newick_str

    def set_name(self, name):
        self.name = name

    def set_len(self, edge_len):
        self.length = float(edge_len)

    def set_parent(self, parent):
        self.parent = parent

    def set_child(self, child):
        if child.parent:
            pass
        self.children.append(child)
        child.parent = self

    def set_sibling(self, sibling):
        self.sibling.append(sibling)

    def set_type(self, cell_type):
        self.cell_type = cell_type

    def inheret(self):
        #self.genome = copy.deepcopy(self.parent.genome)
        #self.alterations = copy.deepcopy(self.parent.alterations)
        self.genome = {}
        chrom_names = list(self.parent.genome.keys())
        for chrom in chrom_names:
            self.genome[chrom] = [[], []]
            for allele in [0, 1]:
                for homolog in self.parent.genome[chrom][allele]:
                    self.genome[chrom][allele].append(homolog.copy())

    def is_leaf(self):
        return len(self.children) == 0
    
    def is_root(self):
        return self.parent is None

    def get_root(self):
        root = self
        while root.parent is not None:
            root = root.parent
        return root

    def add_child(self, name='', edge_len = 0):
        child = Node(name = name, edge_len = float(edge_len), parent = self)
        self.children.append(child)
        return child

    def detach(self):
        if not self.is_root():
            self.parent.children.remove(self)
            self.parent = None

    def iter_postorder(self):
        visit_queue = deque()
        return_queue = deque()
        visit_queue.append(self)

        while visit_queue:
            node = visit_queue.pop()
            return_queue.append(node)
            if not node.is_leaf():
                visit_queue.extend(node.children)
        
        while return_queue:
            node = return_queue.pop()
            yield node

    #level order traversal, i.e. breadth first search from the root. Does include the root itself.
    def iter_descendants(self):
        nodes = deque()
        nodes.append(self)

        while nodes:
            node = nodes.popleft()
            nodes.extend(node.children)
            #if node != self:
            yield node
    
    def iter_preorder(self):
        visit_queue = deque()
        visit_queue.append(self)

        while visit_queue:
            node = visit_queue.pop()
            visit_queue.extend(node.children)
            yield node

    def get_leaves(self):
        return [n for n in self.iter_postorder() if n.is_leaf()]

    def get_height(self):
        if self.is_leaf():
            return 0
        else:
            return max([1 + child.get_height() for child in self.children])

    def get_total_branchlen(self):
        return sum([node.length for node in self.iter_descendants()])

class Tree:
    def __init__(self, root=None, newick=None):
        self.root = root
        self.founder = None
        if newick:
            self.root = str_to_newick(newick)

    def __len__(self):
        return len([leaf for leaf in self.iter_leaves()])

    def print_newick(self, format=0):
        newick_str = self.root.write_newick(format=format)
        print(newick_str)

    def set_founder(self, node):
        self.founder = node
        node.set_type('founder')

    def iter_postorder(self):
        return self.root.iter_postorder()
    
    def iter_descendants(self):
        return self.root.iter_descendants()
    
    def iter_preorder(self):
        return self.root.iter_preorder()
    
    def iter_leaves(self):
        for node in self.iter_descendants():
            if node.is_leaf():
                yield node
    
    def iter_path_to_leaf(self, leaf):
        path = []
        cur_node = leaf
        while not cur_node.is_root():
            path.append(cur_node)
            cur_node = cur_node.parent
        path.append(self.root)
        path.reverse()
        for node in path:
            yield node

    def has_leaf_names(self):
        for leaf in self.iter_leaves():
            if leaf.name == '':
                return False
        return True

    def set_leaf_names(self):
        count = 1
        for leaf in self.iter_leaves():
            if leaf.name == '':
                leaf.name = 'cell' + str(count)
                count += 1

    def set_node_names(self):
        internalcount = 1
        leafcount = 1
        for node in self.iter_descendants():
            if node.is_root():
                if node.cell_type == 'founder':
                    node.name = 'founder'
                else:
                    node.name = 'root'
            elif node.is_leaf():
                #if node.cell_type == 'normal':
                #    node.name = 'normal' + str(leafcount)
                #    leafcount += 1
                #elif node.cell_type == 'pseudonormal':
                #    node.name = 'psd_normal' + str(leafcount)
                #    leafcount += 1
                #else:
                node.name = 'cell' + str(leafcount)
                leafcount += 1
            else:
                if node.cell_type == 'founder':
                    node.name = 'founder'
                else:
                    node.name = 'ancestor' + str(internalcount)
                    internalcount += 1

    def get_tree_height(self):
        return self.root.get_height()

    def get_total_branchlen(self):
        return self.root.get_total_branchlen()

    def ref_node(self, node_name):
        for node in self.iter_descendants():
            if node.name == node_name:
                return node

    def get_mutations(self):
        mutations = []
        for node in self.iter_descendants():
            #mutations += node.whole_chrom_events
            mutations += node.events
        return mutations

    def save(self, file_path, format=0, from_founder=False):
        if from_founder:
            newick_str = self.founder.write_newick(format=format)
        else:
            newick_str = self.root.write_newick(format=format)
        f = open(file_path, 'w+')
        f.write(newick_str)
        f.close()

def str_to_newick(newick_str):
    split_str = newick_str[:-1].split(',')

    cur_node = None
    nodes = []

    for chunk in split_str:
        while chunk[0] == '(':
            new_node = Node()
            if cur_node:
                cur_node.set_child(new_node)
            cur_node = new_node
            chunk = chunk[1:]
        rest = chunk.split(')')
        if ':' in rest[0]:
            idx = rest[0].index(':')
            cur_node.add_child(name=rest[0][:idx], edge_len=rest[0][idx+1:])
        else:
            cur_node.add_child(name=rest[0])
        if len(rest) > 1:
            for part in rest[1:]:
                if ':' in part:
                    idx = part.index(':')
                    cur_node.set_name(part[:idx])
                    cur_node.set_len(part[idx+1:])
                else:
                    cur_node.set_name(part)
                    cur_node.set_len(0)
                if not cur_node.is_root():
                    cur_node = cur_node.parent
    return cur_node

#Tree with random topology by attaching two random leaves. Branch lengths are left undefined.
def gen_tree_random_topology(m, names = None):
    root = Node()
    leaves = [root]

    while len(leaves) < m:
        idx = np.random.randint(0, len(leaves))
        node = leaves.pop(idx)
        c1 = node.add_child()
        c2 = node.add_child()
        leaves.append(c1)
        leaves.append(c2)
    
    count = 0
    for n in root.get_leaves():
        if names:
            n.name = names[count]
        else:
            n.name = 'leaf' + str(count)
        count += 1

    t = Tree(root=root)
    return t

# Creates branch lengths so that each leaf is equally distant from the root.
def add_branchlen_ultrametric(tree):
    node_max_depths = {}
    for node in tree.iter_postorder():
        if node.is_leaf():
            node_max_depths[node] = 1
        else:
            node_max_depths[node] = max([node_max_depths[c] for c in node.children]) + 1
    
    tree_length = float(tree.get_tree_height())
    node_dists = {tree.root: 0.0}
    for node in tree.iter_descendants():
        node.length = (tree_length - node_dists[node.parent]) / node_max_depths[node]
        node_dists[node] = node.length + node_dists[node.parent]
    
# Given a tree topology with existing branch lengths, computes branch lengths uniformly with a random deviation proportional to the tree height.
def add_branchlen_deviation(tree, shape=1):
    for node in tree.iter_descendants():
        x = np.random.gamma(shape)
        node.length = node.length * x

def set_root_branchlen(tree, scale):
    tot = tree.get_total_branchlen()
    tree.root.set_len(tot*scale)

def scale_edge_lengths(tree, place_param):
    leaf_edge_lens = [leaf.length for leaf in tree.founder.iter_leaves()]
    avg_leaf_len = sum(leaf_edge_lens) / len(leaf_edge_lens)
    scalar = place_param / avg_leaf_len

    for node in tree.founder.iter_descendants():
        if node == tree.founder:
            node.length = place_param
        else:
            node.length = node.length * scalar

# OUTDATED
# Calls the ms binary from the given path. 
def call_ms(num_cells, out_path, growth_rate):
    seed_vals = np.random.randint(1, 100000, 3)
    temp_path = os.path.join(out_path, 'temp_log')
    f = open(temp_path, 'w+')
    if growth_rate == 0:
        call = subprocess.call(['mspms', str(num_cells), '1', '-T', '-seeds', str(seed_vals[0]), str(seed_vals[1]), str(seed_vals[2])], stdout=f)
    else:
        call = subprocess.call(['mspms', str(num_cells), '1', '-G', str(growth_rate), '-T', '-seeds', str(seed_vals[0]), str(seed_vals[1]), str(seed_vals[2])], stdout=f)
    f.close()
    with open(temp_path) as f:
        line = f.readline()
        while line[0] != '(':
            line = f.readline()
        tree_str = line[:-1]

    os.remove(temp_path)
    return tree_str

def gen_tree_standard(num_cells, growth_rate):
    demo = ms.Demography()
    demo.add_population(
        initial_size=0.5,
        growth_rate=growth_rate
    )
    ts = ms.sim_ancestry(
        samples=num_cells,
        demography=demo,
        model=ms.StandardCoalescent(),
        ploidy=1
    )
    tree_str = ts.first().as_newick()
    return tree_str

def gen_tree_sweep(num_cells, growth_rate, num_sweep, s):
    demo = ms.Demography()
    demo.add_population(
        initial_size=0.5,
        growth_rate=growth_rate
    )
    sampleset = ms.SampleSet(num_cells, ploidy=1)
    L = 1
    Ne = min(1e8, num_cells*1e4)

    models = []
    for i in range(num_sweep):
        models.append(ms.StandardCoalescent(duration=np.random.uniform()))
        models.append(
            ms.SweepGenicSelection(
                #position=np.random.randint(1, L-1),
                position=0,
                start_frequency=1 / (2*Ne),
                end_frequency=1 - 1 / (2*Ne),
                s = s,
                dt = 1/Ne
            )
        )
    models.append(ms.StandardCoalescent())

    ts = ms.sim_ancestry(
        samples=[sampleset],
        demography=demo,
        model=models,
        sequence_length=L
    )
    tree_str = ts.first().as_newick()
    return tree_str


def make_tumor_tree(tree_type, num_cells, normal_frac, pseudonormal_frac, root_events, out_path, growth_rate, tree_path, num_sweep, s):
    num_normal = round(num_cells * normal_frac)
    if root_events == 0:
        num_pseudonormal = 0
    else:
        num_pseudonormal = round(num_cells * pseudonormal_frac)
    num_aneuploid = num_cells - num_normal - num_pseudonormal

    if num_aneuploid > 0:
        if tree_type == 0:
            if num_sweep > 0:
                tree_str = gen_tree_sweep(num_aneuploid, growth_rate, num_sweep, s)
            else:
                tree_str = gen_tree_standard(num_aneuploid, growth_rate)
            #tree_str = call_ms(num_aneuploid, out_path, growth_rate)
            tree = Tree(newick=tree_str)
        elif tree_type == 1:
            tree = gen_tree_random_topology(num_aneuploid)
        elif tree_type == 2:
            f = open(tree_path)
            tree_str = f.readline()
            f.close()
            tree = Tree(newick=tree_str)
    else:
        tree = Tree()
    
    # Create a new root node and add normal cells as children.
    if num_normal > 0:
        norm_root = Node()
        norm_root.cell_type = 'normal'
        for i in range(num_normal):
            n = Node()
            n.cell_type = 'normal'
            norm_root.set_child(n)
    
    # Handling diverse pseudonormal cells
    if num_pseudonormal > 0:
        # Create a section of the tree for pseudonormals. Distribute them to have anywhere from 1 to 'root_events' number of events. Form the tree accordingly.
        psnorm_root = Node()
        psnorm_root.cell_type = 'pseudonormal'
        placements = Counter(list(np.random.randint(1, root_events+1, size=num_pseudonormal)))
        ps_keys = sorted(list(placements.keys()))

        psnorm_root.events = [None for i in range(ps_keys[0])]
        for i in range(placements[ps_keys[0]]):
            n = Node()
            n.set_type('pseudonormal')
            psnorm_root.set_child(n)

        cur_node = psnorm_root
        cur_mut = ps_keys[0]
        for k in ps_keys[1:]:
            next_int = Node()
            next_int.cell_type = 'pseudonormal'
            next_int.events = [None for i in range(k - cur_mut)]
            cur_mut = k
            for i in range(placements[k]):
                n = Node()
                n.cell_type = 'pseudonormal'
                next_int.set_child(n)
            cur_node.set_child(next_int)
            cur_node = next_int
        if num_aneuploid > 0:
            tree.root.events = [None for i in range(root_events - cur_mut)]

    # Make updates to the tree structure accordingly
    if num_aneuploid > 0:
        tree.set_founder(tree.root)
        if num_normal > 0 and num_pseudonormal == 0:
            norm_root.set_child(tree.root)
            tree.root = norm_root
        elif num_normal == 0 and num_pseudonormal > 0:
            cur_node.set_child(tree.root)
            tree.root = psnorm_root
        elif num_normal > 0 and num_pseudonormal > 0:
            norm_root.set_child(psnorm_root)
            cur_node.set_child(tree.root)
            tree.root = norm_root
        else:
            if root_events == 0:
                tree.root.events = []
            else:
                tree.root.events = [None for i in range(root_events)]
    else:
        if num_normal > 0:
            tree.root = norm_root
            if num_pseudonormal > 0:
                tree.root.set_child(psnorm_root)
        elif num_pseudonormal > 0:
            tree.root = psnorm_root

    return tree

def select_clones(tree, num_clones, criteria, mu, sig):
    ancestral_aneuploids = {}
    for node in tree.founder.iter_descendants():
        if not node == tree.founder and not node.is_leaf():
            ancestral_aneuploids[node] = len(node)

    nodes, sizes = zip(*ancestral_aneuploids.items())

    if criteria == 0:
        if not mu:
            sizes = [s/sum(sizes) for s in sizes]
        else:
            if not sig:
                sig = mu*0.1
            x = scipy.stats.norm(mu, sig)
            sizes = [x.pdf(s) for s in sizes]
            sizes = [s/sum(sizes) for s in sizes]
    elif criteria == 1:
        sizes = [x.length for x in nodes]
        sizes = [x/sum(sizes) for x in sizes]
    
    clone_founders = np.random.choice(nodes, size=num_clones, p=sizes, replace=False)

    for clone in clone_founders:
        clone.cell_type = 'clone'

    return list(clone_founders)
