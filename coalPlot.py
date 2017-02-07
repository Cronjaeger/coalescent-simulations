'''
This is a collection of scripts intended to plot and display coalescents.
'''

import dendropy as dp
from Bio import Phylo
#import networkx
import numpy as np
import finiteSitesModell_investigations as fsmi
from cStringIO import StringIO
import matplotlib.pyplot as plt

class tree_node(object):

    def __init__(self, data = None, label = None, parent = None, children = list(), Tree = None, dist = 1.0):

        if parent is not None and type(parent) is not tree_node:
            raise ValueError('Invalid parent-argument')
        if Tree is not None and type(Tree) is not rooted_tree:
            raise ValueError('Invalid Tree-argument passed.')
        if any(type(c) is not tree_node for c in children):
            raise ValueError('Invalid: Some children passed which were not of type tree_node')

        self.data = data
        self.label = label
        self.Tree = Tree
        self.parent = parent
        self.distanceToParent = dist
        self.children = list(children)

    def __str__(self):
        return str(label)

    def root(self):
        node = self
        while node.parent is not None:
            node = node.parent
        return node

    def leaves(self):
        '''
        returns a list of all leaves subtending the node.
        if the node itself is a leaf, returns [node],
        '''
        output = []
        if len(node.children) == 0:
            output.append(self)
        else:
            for c in self.children:
                output.extend(c.leaves())
        return output

# class tree_edge(object):
#     def __init__(self, origin, terminus, data = None, label = None, length = 1, Tree = None):
#         self.origin = origin
#         self.terminus = terminus
#         self.data = data
#         self.label = label
#         self.length = length
#         self.Tree = Tree

def node_to_string(node):
    return str(node.label)

def node_to_newick(node, label_internal_nodes = False, initialCall = True):
    '''
    prints a newick representation of the phylogeny rooted at 'node'.
    '''

    # we call the function on all chiildren and agregate the results
    out_str = ','.join([node_to_newick(child, label_internal_nodes, False) for child in node.children])

    if len(node.children) > 0:
        out_str = '(%s)'%out_str
        if label_internal_nodes:
            out_str += node_to_string(node)
    else:
        out_str = node_to_string(node)

    if initialCall:
        return '%s;'%out_str
    else:
        return '%s:%f'%(out_str,float(node.distanceToParent))


class rooted_tree(object):

    def __init__(self, root, root_key = 'root'):
        '''
        initialize a rooted tree. The argument root must be of type node.
        '''
        if type(root) is not tree_node:
            raise ValueError('argument \'root\' must have type tree_node')


        root.Tree = self
        self.root = root
        self.nodes = [root]
        self.nodes_dict = {root_key:root}

    def add_node(self, new_node, parent = None, key = None):
        if new_node in self.nodes:
            raise ValueError('Error: atempt to add a node which is already in the tree!')
        if parent is not None:
            if parent.Tree is not self:
                raise ValueError('cannot add node to Tree (specified parent is not in it)')

        new_node.Tree = self

        if parent is not None:
            new_node.parent = parent
            parent.children.append(new_node)

        self.nodes.append(new_node)

        if key == None: ## we were not provided with a key for the node.
            i = len[nodes] #first attempt: the nth node to be added is given the label n
            while self.nodes_dict.has_key(i): # avoid label-collisions.
                i += 1
            #new_node.label = i
            key = i
        self.nodes_dict[key] = new_node

    def subdivide_edge(self, origin, terminus, new_node, d_origin_to_new = 1.0, d_new_to_terminus = 1.0, new_node_key = None):
        '''
        Will insert new_node between origin and terminus
        '''
        try:
            float(d_origin_to_new)
            float(d_new_to_terminus)
        except TypeError:
            raise TypeError('Invalid distances passed.')
        self.add_node(new_node, key = new_node_key)

        #new node replaces terminus in list of children of origin
        i = origin.children.index(terminus)
        origin.children[i] = new_node
        new_node.distanceToParent = d_origin_to_new

        #new node is new parent of terminus
        terminus.parent = new_node
        new_node.children.append(terminus)
        terminus.distanceToParent = d_new_to_terminus


    def str_newick(self,label_internal = True):
        return node_to_newick(self.root,label_internal)

    def str_nexus(self,label_internal = False):
        preamble_str = '#NEXUS\n'
        taxa = [str(node.label) for node in self.nodes]
        taxa_str = '\nBEGIN TAXA;\n    TaxLabels %s;\nEND;\n'%' '.join(taxa)
        tree_str = '\nBEGIN TREES;\n    Tree phylogeny=%s\nEND;\n'%self.str_newick(label_internal)
        return preamble_str + taxa_str + tree_str


def sim_to_tree(sim):
    if type(sim) is not fsmi.simulator_KingmanFiniteSites:
        raise Error('currently only implemented for finite sites coalescent (in other cases, mutations must be handled differently than this method does.)')

    coal = sim.coal

    mutation_events = coal.getMutations()
    mutation_events.sort(reverse = True) #sort with highest time-index first (closest to MRCA)

    merger_events = coal.getChain_with_times()
    merger_events.sort(reverse = True)
    #print merger_events
    # events = coal.getChain_all_events()

    root = merger_events[0]
    root_data = { 't':root[0],
                  'block':root[1].blocks[0],
                  'type':'lineage',
                  'sequence':np.array(sim.MRCA_seq)}
    root_node = tree_node(data = root_data, label='MRCA', dist = None)
    tree = rooted_tree(root_node, root_key = tuple(root[1].blocks[0]))

    # for event in events:
    #     if type(event[1]) is libcoal.partition:
    #         pass
    #     elif type(event[1]) is np.int64:
    #         pass
    #     else
    #         raise Warning('Unknown event-type: %s for event %s'%(str(type(event[1])),str(event)))

    t_old = root[0]
    blocks = [list(root[1].blocks[0])]

    for merger_event in merger_events[1:]:

#        print 'event: t=%f, %s'%(merger_event[0],str(merger_event[1]))

        t_new = merger_event[0]
        blocks_old = list(blocks)

        blocks = list(merger_event[1].blocks)
        blocks_new = [b for b in blocks if b not in blocks_old]
        removed_block = [b for b in blocks_old if b not in blocks][0]

        try:
            assert len(blocks_new) == 2
            assert set(removed_block) == set(blocks_new[0]).union(set(blocks_new[1]))
        except AssertionError:
            print 'blocks:      ',blocks
            print 'blocks_new:  ',blocks_new
            print 'blocks_old:  ',blocks_old
            raise AssertionError()

        parent_node = tree.nodes_dict[tuple(removed_block)]
        t_old = parent_node.data['t']
        assert type(t_old) is not tuple

        block_data_basic = {'t':t_new,'type':'lineage', 'sequence':np.array(parent_node.data['sequence'])}
        for b in blocks_new:
            block_data = dict(block_data_basic)
            block_data['block'] = b
            block_label = label_block(b)
            block_node = tree_node( data = block_data, label = block_label, parent = parent_node, Tree = tree, dist = t_old - t_new)
            tree.add_node(block_node, parent_node, key = tuple(b))

        # we update times in all blocks
        for b in blocks:
            node = tree.nodes_dict[tuple(b)]
            node.data['t'] = t_new
            node.distanceToParent = node.parent.data['t'] - t_new


        for mutation in filter(lambda x: t_new < x[0] and x[0] <= t_old, mutation_events):

            # Extract relevent information from the entry
            t_mutation = mutation[0]
            lineage = mutation[1]
            site = mutation[2]
            shift = mutation[3]
            mutation_index = mutation[4]


            # Determine the parent and child nodes of the lineages in question
            node_below = tree.nodes_dict[tuple(blocks[lineage])]
            node_above = node_below.parent

            # determine the length of the new edges
            t_above = node_above.data['t']
            t_below = node_below.data['t']
            delta_above = t_above - t_mutation
            delta_below = t_mutation - t_below

            #determine how the mutation changes the sequence
            sequence = np.array(node_above.data['sequence'])
            assert all(sequence == node_below.data['sequence'])
            allele_before = sequence[site]
            sequence += shift
            sequence %= 4
            allele_after = sequence[site]

            #update the sequence of the child node
            node_below.data['sequence'] = np.array(sequence)

            node_data = {'t':t_mutation,
                         'type':'mutation',
                         'sequence':sequence,
                         'site':site,
                         'shift':shift,
                         'from_to':(allele_before, allele_after)}

            label = label_mutation(mutation, site, allele_before, allele_after)

            mutation_node = tree_node(data = node_data, label = label)
            mutation_key = 'm%i'%mutation_index

            tree.subdivide_edge( origin = node_above,
                                 terminus = node_below,
                                 new_node = mutation_node,
                                 d_origin_to_new = delta_above,
                                 d_new_to_terminus = delta_below,
                                 new_node_key = mutation_key)

    return tree

def label_mutation(mutation, site, allele_before, allele_after):
    return '%i@%i'%(allele_after,site)
    #return '%ito%i@%i'%(allele_before,allele_after,site)
# def label_mutation(m):
#     return 'x+%i\%4_@site_%i'&(m[4],m[3])

def label_block(b):
    # return '-'.join(map(str,b))
    if len(b) > 1:
        return ''
    else:
        return b[0]

def label_leaf(l):
    return str(l)

def test(n = 10, theta = 1.0, L = 10):
    sim = fsmi.simulator_KingmanFiniteSites(n,theta,L)
    myTree = sim_to_tree(sim)

    # str_newick_test = '(A,(B,C)D);'
    # dpTree_test = dp.Tree.get(data = str_newick_test, schema = 'newick')
    # dpTree_test.print_plot()

    str_newick_sim = myTree.str_newick(True)
    print str_newick_sim
    dpTree_sim = dp.Tree.get(data = str_newick_sim, schema = 'newick')
    dpTree_sim.print_plot()

    phyloTree = Phylo.read(StringIO(str_newick_sim),'newick')
    print phyloTree

    plt.figure()
    Phylo.draw(phyloTree)
    phyloTree.rooted = True
    plt.figure()
    Phylo.draw_graphviz(phyloTree, prog = 'dot')
    plt.draw()

    # str_nexus_sim = tree.str_nexus(True)
    # dpTree_sim_nx = dp.Tree.get(data = str_nexus_sim, schema = 'nexus')
    # dpTree_sim_nx.print_plot()
if __name__ == '__main__':
    #test()
    pass
