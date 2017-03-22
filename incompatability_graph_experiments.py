import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from finiteSitesModell_investigations import *
from TwoSat import *

def incompatability_graph(simulation):
    '''
    take a simulation, compute the incompatibility graph, return a graph-object
    encoding the incompatibility-graph.
    '''

    #Step 1 parse input

    #case: input is a simulated coalescent (used for simulating recurrent mutation)
    if type(simulation) == simulator_KingmanFiniteSites:
        vertices = range(simulation.L)
        edges = simulation.getInconsistentColumnPairs()

    #case: input was generated from parse_ms_output( )
    elif type(simulation) == tuple:
        assert type(simulation[0]) == np.ndarray
        assert type(simulation[1] == list)
        S = simulation[0]
        positions = simulation[1]

        vertices = positions

        incomp_pairs = inconsistentColumnPairs(S)
        edges = [(positions[i],positions[j]) for i,j in incomp_pairs]
    else:
        raise ValueError('unknown simulation-argument')


    #Step 2: generate graph

    # non_trivial_connected_components = fromEdgesToConnectedComponents(edges)
    # trivial_connected_components = [ {v} for v in verticees - set.union(non_trivial_connected_components)]
    G = nx.Graph()
    G.add_nodes_from(vertices)
    G.add_edges_from(edges)

    return G


def simulate_single_graph_and_output_incompatibility_graph(n = 10, L = 20, theta = 20.0, draw_largest_connected_component = False):
    simulation = simulator_KingmanFiniteSites(n,float(theta)/2,L)
    G = incompatability_graph(simulation)
    incompatability_number = len(G.edges())
    highest_degree = max(G.degree().values())
    connected_components_as_sets = [g for g in nx.components.connected.connected_components(G)]
    connected_components_as_graphs = [g for g in nx.components.connected.connected_component_subgraphs(G)]

    largest_connected_component = max(connected_components_as_graphs, key = lambda g: len(g.nodes()))


    # for g_sub in connected_components_as_graphs:
    #     nx.draw_spring(g_sub,with_labels=True, node_color = node_colors(simulation,nodelist = g_sub.nodes()), edge_color = 'black')
    nx.draw_spring(largest_connected_component,with_labels=True, node_color = node_colors(simulation,nodelist = largest_connected_component.nodes()), edge_color = 'black')

    #nx.draw_spring(G,with_labels=True, node_color = node_colors(simulation), edge_color = 'black')
    #plt.show()
    return {    'graph':G,
                'simulation':simulation,
                'components':connected_components_as_graphs,
                'components_as_sets':connected_components_as_sets}

def node_colors(simulation,nodelist = []):
    indexed_palette = {
                    0:'white',
                    1:'pink',
                    2:'turquoise',
                    3:'lightgrey',
                    4:'blue',
                    5:'purple',
                    6:'yellow',
                    7:'darkgreen',
                    8:'lime',
                    9:'teal',
                    10:'black',
                    11:'green',
                    12:'red',
                    13:'orange'}

    if len(nodelist) == 0:
        counts = simulation.site_mutationCount
    else:
        counts = (simulation.site_mutationCount[i] for i in nodelist)

    colors_list = [indexed_palette[min(13,k)] for k in counts]

    return colors_list

def determine_S_of_connected_components(results, print_results = True, return_list = False, with_header = False):
    '''
    Will generate and output the characrters of all columns that contribute to a
    connected component in the incompatibility-graph associated with a simulated
    coalescent with mutation. Each component is treated separately.

    return_list  indicates if the results printed should also be returned (as a list)

    with_header  when set to true will insert a line at the top of each column giving the degree of
                 the node associated with that collumn/character in the binary incompatibility graph.
    '''
    sim = results['simulation']
    S = sim.getS()
    node_lists = [list(s) for s in results['components_as_sets'] if len(s) > 1]
    Graph = results['graph']
    S_list = []

    for column_indices in node_lists:
        column_indices.sort()
        S_reduced = S[:,column_indices]

        if with_header:
            header_f = lambda x: Graph.degree(x) # The header will be the degree associated with a column.
            header = np.array(map(header_f,column_indices),ndmin=2)
            S_reduced = np.r_[header,S_reduced]

        if return_list:
            S_list.append(S_reduced)

        if print_results:
            print S_reduced

    if return_list:
        return S_list

    else:
        return

def run_sim(args):
    p = subprocess.Popen(args,shell = False, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    #print 'run_sim_output: ',output
    return output

def generate_seed():
    '''
    an auxiliary function which allows us to generate a random seed to opass to
    ms, but which depends on the state of the ransom number generator used by
    np.random. This way, seeding the random number generator at the beginning of
    a script with "np.random.set_state()" will also guarantee that the seeds
    used in any sub-processes executring ms will be the same.
    '''
    return int(np.random.sample(1)*(2**31))

def ms_sim_theta(n, theta, N = 1, rho = 0.0 , nsites = 1000000, seed = generate_seed(), with_trees = True, extra_args = []):
    if n<=1 or theta <0 or N<1: raise(ValueError('invalid input to ms_sim_theta: n=%i, theta=%f, N=%i.'%(n,theta,N)))
    if with_trees:
        extra_args.append('-T')
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-t',str(theta),'-rho', str(rho), str(nsites)]+extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)
    return output_parsed

def ms_sim_sites(n, s, N = 1, rho = 0.0, nsites = 1000000, seed = generate_seed(), with_trees = True, extra_args = []):
    if n<=1 or s <1 or N<1: raise(ValueError('invalid input to ms_sim_sites: n=%i, s=%i, N=%i.'%(n,s,N)))
    if with_trees:
        extra_args.append('-T')
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-s',str(s),'-rho', str(rho), str(nsites)] + extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)
    return output_parsed

def parse_ms_output(raw_input):
    'auxiliary function for parsing raw output of ms'
    raw = str(raw_input)
    lines = str.split(raw,'\n')
    chunks =  str.split(raw,r'//')
    # print 'parsing raw input:'
    # print raw
    # print '='*80
    '''
    for refenrence, this is what a chunk looks like:
        ['',
         'segsites: 2',
         'positions: 0.3134 0.5345 ',
         '00',
         '00',
         '01',
         '11',
         '00',
         '',
         '']
    '''
    S_list = []
    trees_list = []
    positions_list = []

    for i,chunk in enumerate(chunks[1:]): #the first entry lists input,seed and trees (if applicable)

        # add an extra empty line to the last line (this way all chunks terminate in two empty lines)
        if chunk == chunks[-1]:
            chunk += '\n'

        #separate the chunk into lines
        chunk_lines = str.split(chunk,'\n')

        #the 0th line is empty

        #the next lines represent trees (if they are simulated)
        if chunk_lines[1][0] == 's':
            n_trees = 0
            trees = []
            #print 'NO trees. chunk_lines[1]=%s'%chunk_lines[1]
        else:
            n_trees = 1
            trees = [parse_tree_from_ms(chunk_lines[n_trees])]
            while chunk_lines[n_trees+1][0] != 's':
                n_trees += 1
                trees.append(parse_tree_from_ms(chunk_lines[n_trees]))

        trees_list.append(trees)

        #the following line contains the number of segregating sites
        segsites_str = chunk_lines[1+n_trees]
        n_segsites = parse_segsites(segsites_str)

        # when there are 0 segregating sites, there is nothing more to parse.
        if n_segsites == 0:
            positions = []
            S = np.array([0],dtype=int,ndmin=2)

        #when there are segregating sites, there is also output to parse
        else:
            positions_str = chunk_lines[2+n_trees]
            positions = parse_positions(positions_str)

            rows_str = chunk_lines[(3+n_trees):-2] # these are the lines which correspond to output
            rows_arr = map(parse_row,rows_str) #each row is converted to a 1d array of integers
            S = np.r_[rows_arr]

        positions_list += [positions]
        S_list.append(S)

    return {'raw':raw,
            'lines':lines,
            'metadata':parse_metadata(chunks[0]),
            'S_list':S_list,
            'positions':positions_list,
            'trees_list':trees_list,
            'n_segsites':n_segsites,
            'experiments':zip(S_list,positions_list,trees_list)}

def parse_segsites(input_str):
    # print 'parsing segsites from string: %s'%input_str
    return int(str.split(input_str)[1])

def parse_positions(input_str):
    # print 'parsing positions from string: %s'%input_str
    return map(float,str.split(input_str)[1:])

def parse_tree_from_ms(input_str):
    '''
    takes a tree-string like '[10](1:0.201,(2:0.078,3:0.078):0.123);', and
    outputs
    '''
    # print 'parsing tree from string: %s'%input_str
    if input_str[0]=='[':
        closing_bracket_position = input_str.index(']')
        sites = int(input_str[1:closing_bracket_position])
        newick_str = input_str[closing_bracket_position+1:]
    elif input_str[0]=='(':
        newick_str = input_str
        sites = 'all'
    else:
        raise TypeError("cannot parse a tree from input string '%s'"%input_str)

    return {'sites':sites,'str':newick_str}

def parse_row(row_str):
    return np.array([int(i) for i in row_str], dtype = int)

def parse_metadata(input_str):
    lines = str.split(input_str, '\n')
    return {'command':lines[0],'seed':int(lines[1])}

def run_experiment_1(N = 1000, sites = 10000, sequences = 10, theta_low_per_site = float(10**-4), theta_high_per_site = float(10**-4), rho_per_site = float(10**-5), single_recombinant_site = True):

    theta_low = sites*theta_low_per_site
    theta_high = sites*theta_high_per_site
    mutation_rate = theta_high/2.0
    rho = rho_per_site * sites

    if single_recombinant_site == True:
        recomb_sites = 2
    else:
        recomb_sites = sites

    # if seed is not None:
    #     np.random.seed(seed)
    # else:
    #     seed = generate_seed()

    #we simulate a coalescent-model with recombination
    recomb_sim = ms_sim_theta(n = sequences, theta = theta_low, N = N, rho = rho, nsites = recomb_sites)
    recomb_experiments = recomb_sim['experiments']


    #we re-set the seed
    #np.random.seed(seed)

    #we simulate with recurrent mutation under the finite sites model
    recurrent_sim = [simulator_KingmanFiniteSites(sequences,mutation_rate,sites,True) for i in xrange(N)]

    #we go over each list of experiments, and generate their incompatibility graphs
    recomb_graphs = map(incompatability_graph,recomb_experiments)
    recurrent_graphs = map(incompatability_graph,recurrent_sim)

    #we find indices of simulations that contain incompatibilities
    incomp_recomb = []
    incomp_recurrent = []
    for i in xrange(N):

        G_recomb = recomb_graphs[i]
        if len(G_recomb.edges()) > 0:
            incomp_recomb.append(i)

        G_recurrent = recurrent_graphs[i]
        if len(G_recurrent.edges()) > 0:
            incomp_recurrent.append(i)

    #generate list of non-trivial incompatibility graphs
    recomb_graphs_incomp = [recomb_graphs[i] for i in incomp_recomb]
    recurrent_graphs_incomp = [recurrent_graphs[i] for i in incomp_recurrent]

    #compute proportion of datasets with incompatabilities
    incomp_rate_recomb = float(len(incomp_recomb))/N
    incomp_rate_recurr = float(len(incomp_recurrent))/N


    if single_recombinant_site:
        #If the following is false, something has gone terribly wrong
        assert all(map(nx.is_bipartite, recomb_graphs_incomp))

    #return a dictionary of computed values for plotting etc.
    return {'recomb_sim':recomb_sim,
            'recurrent_sim':recurrent_sim,
            'recomb_graphs':recomb_graphs,
            'recomb_graphs_incomp':recomb_graphs_incomp,
            'recurrent_graphs':recurrent_graphs,
            'recurrent_graphs_incomp':recurrent_graphs_incomp,
            'incomp_rate_recurr':incomp_rate_recurr,
            'incomp_rate_recomb':incomp_rate_recomb,
            'recomb_index_incomp':incomp_recomb,
            'recurr_index_incomp':incomp_recurrent
            }

def summary(S):
    '''
    Computes a range of summary-statistics for a dataset encoded in an S-matrix
    '''
    n_sequences = S.shape[0]
    n_chars = S.shape[1]

    partitions = map(to_partition,np.transpose(S))

    state_distribution = {1:0,2:0,3:0,4:0}
    n_blocks_in_char = map(len, partitions)
    for k in n_blocks_in_char:
        state_distribution[k] += 1
    assert sum(state_distribution.values()) == n_chars

    block_sizes = map(lambda x: map(len,x), partitions)
    block_sizes_sorted = []
    for integer_composition in block_sizes:
        integer_partiton = list(integer_composition)
        integer_partiton.sort(reverse = True)
        block_sizes_sorted.append(integer_partiton)

    informative_columns = get_informative_collumns(S)
    n_informative_columns = len(informative_columns)

    G = incompatability_graph((S,range(S.shape[1])))
    stariness = compute_staryness(G)
    ranked_handshake_dist = ranked_handshake_distribution(G)

    incompatible_groups = reduce(lambda x,y: x+y, [list(s) for s in nx.connected_components(G) if len(s) > 1], [])

    S_only_incomplatible = S[:,incompatible_groups]

    summary_statistics = \
    {
    'data':np.array(S),
    'n seq':n_sequences,
    'n chars':n_chars,
    'pairwise incompatibility graph':G,
    'stariness':stariness,
    'ranked_handshake_dist':ranked_handshake_dist,
    'chars as partitions on taxa':partitions,
    'state distribution':state_distribution,
    'n informative chars':n_informative_columns,
    'S with only problematic characters':S_only_incomplatible
    }


    return summary_statistics

def columns_incompatible_with_tree(simulation):
    '''returns the list of indicees for all columns in the S-matrix of simulation
which are incompatible with the underlying tree'''
    #assert isinstance(simulation,simulator_KingmanFiniteSites)

    partition_chain = simulation.coal.getPartitionChain()

    # go from [ [.], [.], ... , [.] ] to [ set([.]), set([.]), ... , set([.]) ]
    partition_chain = map(lambda partition: map(set, partition.blocks), partition_chain)

    S = simulation.getS()
    #informative_columns = get_informative_collumns(S)
    #chars_as_partitions = np.apply_along_axis(func1d = to_partition, axis = 1, arr = S)
    chars_as_partitions = [ to_partition( S[:,i] ) for i in range(S.shape[1])]
    informative_indices = [i for i in range(S.shape[1]) if is_informative(S[:,i])]

    compatible_informative_indices = list(informative_indices)
    for partition in partition_chain:
        compatible_informative_indices = filter(lambda i: partition_compatibility(partition, chars_as_partitions[i]), compatible_informative_indices)

    incompatible_indices = [i for i in informative_indices if i not in compatible_informative_indices]

    return incompatible_indices

def block_compatibility(block_1,block_2):
    assert type(block_1) == type(block_2) == set
    return block_1.issubset(block_2) or block_1.issuperset(block_2) or block_1.isdisjoint(block_2)

def partition_compatibility(part1,part2):
    return all([all(map(lambda B: block_compatibility(A,B), part2)) for A in part1])

def generate_plots_experiment_1(results, rows_max = 5, savepath = './', extensions = ['.png','.jpeg','.svg','.pdf'], base_node_size = 20,alpha = 0.2, drop_deg_0_nodes = True, drop_non_segregating_sites = True):

    #for politting
    cols = 2
    rows = rows_max

    recurrent_graphs_incomp = results['recurrent_graphs_incomp']
    recomb_graphs_incomp = results['recomb_graphs_incomp']


    recurrent_graphs_to_plot = recurrent_graphs_incomp[:(min(rows_max,len(recurrent_graphs_incomp)))]
    recomb_graphs_to_plot = recomb_graphs_incomp[:(min(rows_max,len(recomb_graphs_incomp)))]

    rows = max(len(recurrent_graphs_to_plot),len(recomb_graphs_to_plot))

    fig = plt.figure(figsize = (16,5*rows))

    for i,G in enumerate(recurrent_graphs_to_plot):
        G = nx.Graph(G) # copy the graph so we don't change anything in results

        #determine auxiliary statistics
        sim_index = results['recurr_index_incomp'][i]
        d_summary = summary(results['recurrent_sim'][sim_index].getS())

        large_components = filter(lambda x: len(x)>1, nx.connected_components(G))
        large_component_sizes = map(len,large_components)
        large_component_sizes_str = ' '.join(map(str,large_component_sizes))

        #sumarize aux. statistics
        title_str = 'chars: %i (%i inf.), state-dist: %s\nsize: %s,$\\star$ = %0.2f'%(d_summary['n chars'],d_summary['n informative chars'],str(tuple(d_summary['state distribution'].values())),large_component_sizes_str,d_summary['stariness'])

        sites = results['recurrent_sim'][i].L


        char_states_at_site = map(lambda k: len(set(results['recurrent_sim'][sim_index].sequences[:,k])),range(sites))
        palette = {1:'black',2:'blue',3:'red',4:'green'}
        color_dict = dict([(node,palette[char_states_at_site[node]]) for node in G.nodes()])

        if drop_deg_0_nodes:
            remove_nodes_with_degree_0(G)
        elif drop_non_segregating_sites:
            irrelevant_sites = filter(lambda k: char_states_at_site[k] == 1, range(sites))
            #print len(irrelevant_sites)
            G.remove_nodes_from(irrelevant_sites)


        plt.subplot(rows,cols,2*i+1)
        plt.title(title_str)
        # nx.draw_spring(G, node_color = 'red')
        draw_FSM_graph(G, base_node_size = base_node_size,alpha = alpha, color_dict = color_dict)

    for i,G in enumerate(recomb_graphs_to_plot):

        G = nx.Graph(G) # copy the graph so we don't change anything in results
        if drop_deg_0_nodes:
            remove_nodes_with_degree_0(G)

        large_component_graphs = filter(lambda x: len(x)>1, nx.connected_component_subgraphs(G))
        size_strs = []
        for g in large_component_graphs:
            if nx.is_bipartite(g):
                l,r = nx.bipartite.sets(g)
                size_strs.append('%i-%i'%(len(l),len(r)))
        size_str = ' '.join(size_strs)
        #size_l,size_r = sizes[0],sizes[1]

        sim_index = results['recomb_index_incomp'][i]
        d_summary = summary(results['recomb_sim']['S_list'][sim_index])
        title_str = 'chars: %i (%i inf.), blocks: %s, $\\star$ = %0.2f'%(d_summary['n chars'],d_summary['n informative chars'],size_str,d_summary['stariness'])


        plt.subplot(rows,cols,2*(i+1))
        plt.title(title_str)
        # nx.draw_spring(G, node_color = 'blue')
        draw_IS_w_recomb_graph(G, base_node_size = base_node_size,alpha = alpha)

    # graphs = rows*cols
    # for i in range(graphs):
    #     plt.subplot(rows,cols,i)
    #     #print i
    #     results = simulate_single_graph_and_output_incompatibility_graph(n =3, L = 100, theta = 10.0, draw_largest_connected_component = True)
    #     #print
    #     determine_S_of_connected_components(results,print_results=False,with_header=True)
    plt.draw()

    return fig
    # sample  incompatible datasets from each experiment, and plot them

def draw_FSM_graph(G, alpha = .2, base_node_size = 10, color_dict = None):

    sizeDict = size_proportional_to_degree(G, base = base_node_size)
    nodelist = sizeDict.keys()
    sizes = sizeDict.values()
    if color_dict is None:
        node_color = 'red'
    else:
        node_color = [color_dict[node] for node in nodelist]
    #nx.draw_spring(G, node_color = 'red')
    #pos = nx.spring_layout(G)
    pos = node_positions(G)
    nx.draw(G, nodelist = nodelist, width = 0.4, pos = pos, node_color = node_color, alpha = alpha, node_size = sizes, linewidths = 0.0)

def draw_IS_w_recomb_graph(G, alpha = .2, base_node_size = 10, labels = {},font_size = 12):

    sizeDict = size_proportional_to_degree(G, base = base_node_size)
    nodelist = sizeDict.keys()
    sizes = sizeDict.values()

    if nx.bipartite.is_bipartite(G):
        #pos = bipartite_layout(G)
        posMap = bipartite_layout
    else:
        #pos = nx.spring_layout(G)
        posMap = nx.spring_layout
    pos = node_positions(G,componentMap = posMap)
    nx.draw(G, pos = pos, node_color = 'blue', width = 0.2, nodelist = nodelist, node_size = sizes, alpha = alpha, linewidths = 0.0, labels = labels, font_size = font_size)
    #nx.draw(G, pos = pos, node_color = G.nodes(),vmin = 0.5,vmax = 0.6, width = 0.2, nodelist = nodelist, node_size = sizes, alpha = alpha)


def node_positions(G,componentMap = nx.spring_layout, whitespace = 0.2):
    assert type(G) is nx.Graph
    assert 0 <= whitespace < 1.0

    subgraphs = [g for g in nx.connected_component_subgraphs(G)]
    n_components = len(subgraphs)

    if n_components == 1:
        return componentMap(G)

    #widths = [1.0/len(subgraphs) for i in range(len(subgraphs))]
    widths = [(1 - whitespace)*float(len(g))/len(G) for g in subgraphs]
    #print sum(widths),widths

    positions = {}

    for i,g in enumerate(subgraphs):
        pos_component = componentMap(g)
        shift = sum(widths[:i]) + i*whitespace/n_components
        for node in pos_component:
            x,y = pos_component[node]
            positions[node] = (shift + (widths[i]*x), y)

    return positions


# def color_prop_to_pos(G):
#     map(float,G.nodes())

def size_proportional_to_degree(G,base = 10,deg_0_size = 5):
    #hd = handshake_distribution(G):
    # sizes = np.zeros(len(G),dtype = int)
    #
    # for i,deg in enumerate(G.degree().values()):
    #     sizes[i] = base*deg
    #
    # return sizes
    sizeDict = G.degree()
    for v in sizeDict:
        sizeDict[v] += 1
        sizeDict[v] *= base
    return sizeDict


def bipartite_layout(G):

    assert nx.bipartite.is_bipartite(G)

    left,right = nx.bipartite.sets(G)

    size_left = len(left)
    size_right = len(right)

    positions = {}

    for i,node in enumerate(left):
        positions[node] = (float(0),float((size_left - 2*i))/size_left)

    for i,node in enumerate(right):
        positions[node] = (float(1),float((size_right - 2*i))/size_right)
    # for i,node in enumerate(left):
    #     positions[node] = (float(-1),float((size_left - 2*i))/size_left)
    #
    # for i,node in enumerate(right):
    #     positions[node] = (float(1),float((size_right - 2*i))/size_right)

    return positions

# def generate_all_k_tuples_of_i_state_characters_on_n_sequences(k,i,n):
#     out = []
#     states = range(i)
#     for

# def generate_plots_experiment_2(results):
#     '''
#     Will analyze the degree-distribution
#     '''
#     #compute the average degree_distribution
#
#     #compute the averate degree_distribution conditional on incompatibility
#
#     #compute the average value of the degree of the most central node.
#     pass

def compute_staryness(G):
    '''
    Takes a graph as input, and returns the quantity deg_max/|E|
    where deg_max is the maximum degree of a node in the graph, and
    |E| is the total number of edges in the graph.
    '''
    if len(G.edges()) == 0:
        return 0
    return float(max(G.degree().values()))/len(G.edges())

def handshake_distribution(G):
    '''
    returns a dictionary with entries (v,p(v)) where p(v)=deg(v)/2|E|, v is a
vertex of G and |E| is the number of edges of V
    '''
    d = G.degree() # a dictionary of the form (node, degree)
    deg_sum = sum(d.values()) # equals twice the number of edges by the handshake-lemma
    if deg_sum == 0:
        for key in d.keys():
            d[key] = 0.0
    else:
        deg_sum = float(deg_sum)
        for key in d.keys():
            d[key] /= deg_sum
    return d

def ranked_handshake_distribution(G):
    '''
    returns the values of the handshake-distribution sorted in non-ascending order
    '''
    d = handshake_distribution(G)
    s = d.values()
    s.sort(reverse=True)
    return s


def remove_nodes_with_degree_0(G):
    assert type(G) == nx.Graph
    G.remove_nodes_from([x for x,y in G.degree().items() if y == 0])

def generate_random_characters(n_chars,n_states,n_sequences):
    S_transposed = np.random.randint(n_states, size = (n_chars,n_sequences))
    return map(list,S_transposed)

def general_compatibility_test(chars):
    if len(chars) == 2:
        return two_char_compatibility_test(chars)
    elif max(map(lambda x: len(set(x)), chars)) <= 3:
        return three_char_compatibility_test(chars)
    else:
        return NotImplemented

        # G = nw.Graph(partition_intersection_graph(chars)['E'])
        # return restricted_chordal_completion_exists(G)

def minimal_inconsistent_tripples(S, ancestral_type_known = True, inconsistent_pairs = None, disregard_non_informative = False):
    if inconsistent_pairs is None:
        inconsistent_pairs = inconsistentColumnPairs(S,ancestral_type_known)
    inconsistent_pairs = set(inconsistent_pairs) #set supports faster look-up than list

    min_triples = []

    if disregard_non_informative:
        informative = filter(lambda i: is_informative(S[:,i],ancestral_type_known), range(S.shape[1]))
    else:
        informative = range(S.shape[1])
    #m = len(informative)
    for i,s1 in enumerate(informative):
        for j,s2 in enumerate(informative[(i+1):],start = i+1):
            if (s1,s2) not in inconsistent_pairs:
                for s3 in informative[(j+1):]:
                    if (s1,s3) not in inconsistent_pairs and (s2,s3) not in inconsistent_pairs:
                        if three_state_compatability_test((S[:,s1],S[:,s2],S[:,s3])):
                            min_triples.append((s1,s2,s3))
    return min_triples

def three_state_compatability_test(chars):
    '''Implementing the compatibility test for 3-state characters of Dress and Steel
1992'''

    # Each char with three states is turned into a list of three chars.
    # 1 & 2 state chars are turned into a lsit with that same char as its only
    # element

    #triples = map(three_state_to_binary_tripple, chars)
    triples = map(char_to_binary_tripple, chars) # we turn every character into 3 binary characters

    I = range(len(triples))
    #assert len(I) == (len(chars)) #assert that every character has been extended

    #we deine an S-matrix of extended characters
    n_chars_extended = sum(map(len,triples))
    n_seq = len(chars[0])
    S_extended = np.zeros((n_seq,n_chars_extended),dtype = int)

    I_to_J = dict([(i,[]) for i in I])
    J_to_I = dict()

    #we fill out the entries in that S-matrix
    j = 0
    for i,extended_chars in enumerate(triples):
        for extended_char in extended_chars:
            S_extended[:,j] = extended_char
            I_to_J[i].append(j)
            J_to_I[j] = i
            j += 1

    incompatible_pairs = inconsistentColumnPairs(S_extended)

    if len(incompatible_pairs) == 0:
        #the entire set of extended chars is compatible
        return True
    else:
        # J_to_J = dict([(j,[k for k in I_to_J[J_to_I[j]] if k != j]) for j in J_to_I.keys()])

        K = find_consistent_subset(incompatible_pairs = incompatible_pairs, constraint_sets = I_to_J.values())

        if K is None:
            return False
        else:
            return True

def find_consistent_subset(incompatible_pairs,constraint_sets):
    '''Tries to find a subset of the columns of S with the constraint that at
least one element from each incompatible pair is not included, but at most one
element of each constaraint-set is excluded (a constraint-set of smaller size
must have all entries included)'''
    # if len(incompatible_pairs) == 0:
    #     pass #TODO: finish this!
    #     #return constraint_dict
    # else:
    #     #excluded = {}
    #     for a,b in incompatible_pairs:
    #         pass

    #verify that all constraint-sets have length 3 or 1
    assert all(map(lambda length: length == 1 or length == 3, map(len,constraint_sets)))

    #these represent characters with only two states. they MUST be included.
    non_extended_chars = [list(x)[0] for x in constraint_sets if len(x) == 1]

    #In the current implementation all chars should be extended. we verify this.
    assert len(non_extended_chars) == 0

    # these represent groups of three binary characters created from just one
    # character
    extended_chars = [list(x) for x in constraint_sets if len(x) == 3]

    # pruned_extended_chars = list(extended_chars)
    # for i in non_extended_chars:
    #     pruned_extended_chars = [filter(lambda j: j!= i,chars) for chars in pruned_extended_chars]
    #
    # if min(map(len,len_3_elements)) < 2:
    #     # no consistent subsdet exists, since 2 out of 3 chars in each
    #     # trippe must be included.
    #     return []

    #we first construct all clauses that comprize the 2-SAT problem
    clauses_from_incompatible_pairs = inconsistent_pairs_to_bool_var_pairs(incompatible_pairs)
    clauses_from_extended_chars =  flatten(map(constraint_to_bool_var_pairs,extended_chars))
    clauses_all = clauses_from_incompatible_pairs + clauses_from_extended_chars

    #we then solve the associated two-sat-problem
    twoSAT_problem = TwoSatSolver(clauses_all)
    K = twoSAT_problem.solution_trues()
    K.sort()

    #chars_to_include = non_extended_chars + K
    #chars_to_include.sort()

    return K

def flatten(LOL):
    '''Will lake LOL (a List Of Lists)
    = [[x_1, x_2, ... , x_k], [y_1, y_2, ..., y_r], ... ]
    and return
    [x_1,x_2, ... , x_k, y_1, y_2, ... , y_r, ...]'''
    return reduce(lambda x,y: x+y, LOL, [])

def inconsistent_pairs_to_bool_var_pairs(pair_list):
    #return reduce(lambda x,y: x+[(bool_var(y[0],True),bool_var(y[1],False)), (bool_var(y[0],False),bool_var(y[1],True))], pair_list, [])
    return [(bool_var(x,False),bool_var(y,False)) for x,y in pair_list]

def char_to_binary_tripple(char):
    '''Takes any character with at most 3 states (encoded  as a vector), and
    returns a list of two state characters (used in the reduction to 2-SAT)
    outlined in chapter 14 of Gusfeld 2003.
    A two or one state character will be padded with characters of the form
    [0,0,...0], to yield 3 characters total.
    '''
    states = set(char)
    if len(states) <= 3:
        return [ [int(y == x) for y in char] for x in states] + [ len(char)*[0] for i in range(3-len(states))]
    else:
        raise TypeError('char has more than 3 states!')


def three_state_to_binary_tripple(char):
    states = set(char)
    if len(set(char)) < 3:
        return [char]
    elif len(states) == 3:
        states = set(char)
        return [ [int(y == x) for y in char] for x in states]
    else:
        raise TypeError('char has more than 3 states!')

def constraint_to_bool_var_pairs(constraint):
    '''
    takes a collection of at most three indices and returns:
    [(i ,j), (i, k), (j, k)] if constraint = (i,j,k)
    [(i,j)] if constraint (i,j)
    [(i,i)] if constraint (i,i)
    '''
    constraint = map(lambda i: bool_var(i,True), constraint)
    if len(constraint) > 3:
        raise ValueError('input must be of length at most 3')
    elif len(constraint) == 3:
        x = constraint
        return [ (x[0], x[1]),  (x[0], x[2]),  (x[1], x[2]) ]

    elif len(constraint) == 2:
        return [ (x[0], x[1]) ]
    elif len(constraint) == 1:
        return[(x[0],x[0])]
    else:
        raise ValueError('input must be of length at least 1')



###
# The 2-SAT solver has been mooved to a separate file.
###
# class TwoSatSolver(object):
#
#     def __init__(self,clauses):
#         '''Solve 2-SAT given clauses in conjunctive normal form
# Clauses should be a list of pairs of bool_var-instances. Each pair (A, B) is
# interpreted as the statement "A or B".
#  Given input [(A1, B1), (A2,B2), ... ] we
# solve the associated 2-SAT problem which in conjunctive normal form is
# expressed: (A1 or B1) and (A2 or B2) and (A3 or B3) ... '''
#         if any(map(lambda clause: len(clause) != 2, clauses)):
#             raise ValueError('clauses must be a list of pairs')
#         if not all(map(lambda clause: isinstance(clause[0],bool_var) and isinstance(clause[1],bool_var), clauses)):
#             raise TypeError('All clauses must be pairs of instances of bool_var')
#
#         self.clauses = list(clauses)
#
#         self.implications = flatten(map(clause_to_edge,self.clauses))
#
#         '''D is an implication graph. each implication of the form A => B is
#         represented as a directed edge from A to B'''
#         self.D = nx.DiGraph()
#         self.D.add_edges_from(self.implications)
#
#         self.solvable = None
#         self.solution = []
#
#         self.scc = nx.components.strongly_connected_components(self.D)
#         for c in self.scc:
#             c_sorted = list(c)
#             c_sorted.sort()
#             ''' variables are now sorted by label; It is therefore sufficient to
#             verify that c_sorted[i] != c_sorted[i+1] holds for all
#             '''
#             for i in range(len(c) - 1):
#                 if c_sorted[i] == c_sorted[i+1]:
#                     self.sollvable = False
#
#         solution_indices = []
#         if self.solvable is None:
#             self.solvable = True
#             C = nx.condensation(self.D)
#
#             for v in nx.topological_sort(C):
#                 literals = C.node[v]['members']
#                 if set([x.index() for x in literals]).isdisjoint(set(solution_indices)):
#                     for x in literals:
#                         self.solution.append(x)
#                         solution_indices.append(x.index())
#
#             self.solution.sort()
#             assert set(x.index() for x in self.solution) == set([x.index() for x in self.D.nodes()])
#
#     def asDict:
#         return dict([(x.name,x.value) for x in self.solution])
#
#     def solution_trues:
#         return [x.name for x in self.solution if x.value is True]
#
#     def solution_false:
#         return [x.name for x in self.solution if x.value is False]
#
# def flatten(LOL):
#     '''Will lake LOL (a List Of Lists)
#     = [[x_1, x_2, ... , x_k], [y_1, y_2, ..., y_r], ... ]
#  and return
#     [x_1,x_2, ... , x_k, y_1, y_2, ... , y_r, ...]'''
#     return reduce(lambda x,y: x+y, LOL, [])
#
# class bool_var(object):
#     '''A class to represent boolean variables for use when solving 2-SAT.
# bool_var(i,True) = "x_i"; bool_var(i,False) = "not x_i";
# i should be a hashable type
#     '''
#     def __init__(self,name,value = True):
#         #assert type(value) is bool
#         self.name = name
#         self.value = bool(value)
#
#     def __hash__(self):
#         return hash((self.name,self.value))
#
#     def __str__(self):
#         return 'x_%s'%str(self.name) if self.value else '-x_%s'%str(self.name)
#
#
#     def __repr__(self):
#         #return 'boolean variable(%s)'%str(self)
#         return str(self)
#
#     def __eq__(self,other):
#         if not isinstance(other, bool_var):
#             return False
#         else:
#             return self.name == other.name and self.value == other.value
#
#     def __cmp__(self,other):
#         if not isinstance(other,bool_var):
#             return cmp(self.name,other)
#         else:
#             return cmp(self.name,other.name)
#
#     def negated(self):
#         return bool_var(self.name, not self.value)
#
#     def getName(self):
#         return self.name
#
#     def index(self):
#         '''alias for getName'''
#         return self.getName()
#
#     def as_bool(self):
#         return self.value
#
# def neg(bool_var):
#     return bool_var.negated()
#
# def inconsistent_pairs_to_bool_var_pairs(pair_list):
#     #return reduce(lambda x,y: x+[(bool_var(y[0],True),bool_var(y[1],False)), (bool_var(y[0],False),bool_var(y[1],True))], pair_list, [])
#     return [(bool_var(x,False),bool_var(y,False)) for x,y in pair_list]
#
# def constraint_to_bool_var_pairs(constraint):
#     '''
#     takes a collection of at most three indices and returns:
#     [(i ,j), (i, k), (j, k)] if constraint = (i,j,k)
#     [(i,j)] if constraint (i,j)
#     [(i,i)] if constraint (i,i)
#     '''
#     constraint = map(lambda i: bool_var(i,True), constraint)
#     if len(constraint) > 3:
#         raise ValueError('input must be of length at most 3')
#     elif len(constraint) == 3:
#         x = constraint
#         return [ (x[0], x[1]),  (x[0], x[2]),  (x[1], x[2]) ]
#
#     elif len(constraint) == 2:
#         return [ (x[0], x[1]) ]
#     elif len(constraint) == 1:
#         return[(x[0],x[0])]
#     else:
#         raise ValueError('input must be of length at least 1')
#
# def clause_to_edge(clause):
#     assert len(clause) == 2
#     assert type(clause[0]) == type(clause[1]) == bool_var
#     x,y = clause
#     return [(neg(x),y) , (neg(y),x)]
#
# def char_to_binary_tripple(char):
#     '''Takes any character with at most 3 states (encoded  as a vector), and
#     returns a list of two state characters (used in the reduction to 2-SAT)
#     outlined in chapter 14 of Gusfeld 2003.
#     A two or one state character will be padded with characters of the form
#     [0,0,...0], to yield 3 characters total.
#     '''
#     states = set(char)
#     if len(states) <= 3:
#         return [ [int(y == x) for y in char] for x in states] + [ len(char)*[0] for i in range(3-len(states))]
#     else:
#         raise TypeError('char has more than 3 states!')
#
#
# def three_state_to_binary_tripple(char):
#     states = set(char)
#     if len(set(char)) < 3:
#         return [char]
#     elif len(states) == 3:
#         states = set(char)
#         return [ [int(y == x) for y in char] for x in states]
#     else:
#         raise TypeError('char has more than 3 states!')

def restricted_chordal_completion_exists(G):
    '''
    NOT IMPLEMENTED!!!
    '''
    assert type(G) == nw.Graph

    if nw.chordal.is_chordal(G):
        return True
    else:
        return NotImplemented
        #  C = long_cycles_iterator(G)
        #  restricted_chordalization_found = False
        #  while C.next() is not None:
        #     pass
            #return restricted_chordal_completion_exists()

# def ms_sim(arg_string):
#     command = 'ms' + arg_string
#     p = subprocess.Popen(command,Shell = True, stdout=subprocess.PIPE)
#     output = p.communicate()[0]
#     #p.kill()
#     return output
