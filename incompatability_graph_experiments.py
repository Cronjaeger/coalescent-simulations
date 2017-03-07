import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from finiteSitesModell_investigations import *

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
        title_str = 'chars: %i (%i inf.), state-dist: %s, size: %s,$\\star$ = %0.2f'%(d_summary['n chars'],d_summary['n informative chars'],str(tuple(d_summary['state distribution'].values())),large_component_sizes_str,d_summary['stariness'])

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

    G = nw.Graph(partition_intersection_graph(chars)['E'])

    return restricted_chordal_completion_exists(G)

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
