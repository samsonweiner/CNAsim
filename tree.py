import numpy as np
import subprocess
from collections import deque
import copy

class Node:
    def __init__(self, name='', edge_len = 0, parent = None):
        self.name = name
        self.length = edge_len
        self.parent = parent
        self.children = []
        self.events = []
        #self.whole_chrom_events = []
        self.genome = None
        #self.alterations = None
        self.profile = None

    def __str__(self):
        if self.name:
            return str(self.name)
        else:
            return ''
    
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
        self.children.append(child)
        child.parent = self

    def set_sibling(self, sibling):
        self.sibling.append(sibling)

    def inheret(self):
        self.genome = copy.deepcopy(self.parent.genome)
        #self.alterations = copy.deepcopy(self.parent.alterations)

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
        if newick:
            self.root = str_to_newick(newick)
        self.cell_names = [node.name for node in self.iter_leaves()]

    def print_newick(self, format=0):
        newick_str = self.root.write_newick(format=format)
        print(newick_str)

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
                leaf.name = 'leaf' + str(count)
                count += 1
        self.cell_names = [node.name for node in self.iter_leaves()]

    def set_node_names(self):
        internalcount = 1
        leafcount = 1
        for node in self.iter_descendants():
            if node.is_root():
                node.name = 'root'
            elif not node.is_leaf():
                node.name = 'internal' + str(internalcount)
                internalcount += 1
            elif node.is_leaf():
                node.name = 'leaf' + str(leafcount)
                leafcount += 1
        self.cell_names = [node.name for node in self.iter_leaves()]

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

    def save(self, file_path, format=0):
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
def gen_random_topology(m, names = None):
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

# Calls the ms binary from the given path.
def call_ms(ms_path, num_cells, out_path, growth_rate):
    seed_vals = np.random.randint(1, 100000, 3)
    f = open(out_path + 'temp_log', 'w+')
    if growth_rate == 0:
        call = subprocess.call([ms_path, str(num_cells), '1', '-T', '-seeds', str(seed_vals[0]), str(seed_vals[1]), str(seed_vals[2])], stdout=f)
    else:
        call = subprocess.call([ms_path, str(num_cells), '1', '-G', str(growth_rate), '-T', '-seeds', str(seed_vals[0]), str(seed_vals[1]), str(seed_vals[2])], stdout=f)
    f.close()
    with open(out_path + 'temp_log') as f:
        line = f.readline()
        while line[0] != '(':
            line = f.readline()
        tree_str = line[:-1]

    call = subprocess.call(['rm', out_path + 'temp_log'])
    #f = open(out_path + 'tree.nwk', 'w+')
    #f.write(tree_str)
    #f.close()
    return tree_str