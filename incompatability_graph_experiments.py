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

def ms_sim_theta(n, theta, N = 1, rho = 0.0 , nsites = 1000000, seed = generate_seed(), extra_args = []):
    if n<=1 or theta <0 or N<1: raise(ValueError('invalid input to ms_sim_theta: n=%i, theta=%f, N=%i.'%(n,theta,N)))
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-t',str(theta),'-rho', str(rho), str(nsites)]+extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)
    return output_parsed

def ms_sim_sites(n, s, N = 1, rho = 0.0, nsites = 1000000, seed = generate_seed(), extra_args = []):
    if n<=1 or s <1 or N<1: raise(ValueError('invalid input to ms_sim_sites: n=%i, s=%i, N=%i.'%(n,s,N)))
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-s',str(s),'-rho', str(rho), str(nsites)] + extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)
    return output_parsed

def parse_ms_output(raw_input):
    'auxiliary function for parsing raw output of ms'
    raw = str(raw_input)
    lines = str.split(raw,'\n')
    chunks =  str.split(raw,r'//')
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
    positions_list = []

    for chunk in chunks[1:]: #the first entry lists input,seed and trees (if applicable)

        # add an extra empty line to the last line (this way all chunks terminate in two empty lines)
        if chunk == chunks[-1]:
            chunk += '\n'

        #separate the chunk into lines
        chunk_lines = str.split(chunk,'\n')

        #the 0th line is empty

        #the 1st line contains the number of segregating sites
        segsites_str = chunk_lines[1]
        n_segsites = parse_segsites(segsites_str)

        # when there are 0 segregating sites, there is nothing more to parse.
        if n_segsites == 0:
            positions = []
            S = np.array([0],dtype=int,ndmin=2)

        #when there are segregating sites, there is also output to parse
        else:
            positions_str = chunk_lines[2]
            positions = parse_positions(positions_str)

            rows_str = chunk_lines[3:-2] # these are the lines which correspond to output
            rows_arr = map(parse_row,rows_str) #each row is converted to a 1d array of integers
            S = np.r_[rows_arr]

        positions_list += [positions]
        S_list.append(S)

    return {'raw':raw,
            'lines':lines,
            'metadata':parse_metadata(chunks[0]),
            'S_list':S_list,
            'positions':positions_list,
            'n_segsites':n_segsites,
            'experiments':zip(S_list,positions_list)}

def parse_segsites(input_str):
    return int(str.split(input_str)[1])

def parse_positions(input_str):
    return map(float,str.split(input_str)[1:])

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

    if single_recombinant_site:
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

    #compute proportion of datasets with incompatabilities
    incomp_rate_recomb = float(len(incomp_recomb))/N
    incomp_rate_recurr = float(len(incomp_recurrent))/N

    if single_recombinant_site:
        #to check for wierdness
        assert all(map(nx.is_bipartite, recomb_graphs))

    #return a dictionary of computed values for plotting etc.
    return {'recomb_sim':recomb_sim,
            'recurrent_sim':recurrent_sim,
            'recomb_graphs':recomb_graphs,
            'recomb_graphs_incomp':[recomb_graphs[i] for i in incomp_recomb],
            'recurrent_graphs':recurrent_graphs,
            'recurrent_graphs_incomp':[recurrent_graphs[i] for i in incomp_recurrent],
            'incomp_rate_recurr':incomp_rate_recurr,
            'incomp_rate_recomb':incomp_rate_recomb
            }

def generate_plots_experiment_1(results, rows_max = 5, savepath = './', extensions = ['.png','.jpeg','.svg','.pdf']):
    cols = 2
    rows = rows_max

    recurrent_graphs_incomp = results['recurrent_graphs_incomp']
    recomb_graphs_incomp = results['recomb_graphs_incomp']

    recurrent_graphs_to_plot = recurrent_graphs_incomp[:(min(5,len(recurrent_graphs_incomp)))]
    recomb_graphs_to_plot = recomb_graphs_incomp[:(min(5,len(recomb_graphs_incomp)))]

    rows = max(len(recurrent_graphs_to_plot),len(recomb_graphs_to_plot))

    fig = plt.figure(figsize = (10,4*rows))

    for i,G in enumerate(recurrent_graphs_to_plot):

        G = nx.Graph(G) # copy the graph so we don't change anything in results
        remove_nodes_with_degree_0(G)

        plt.subplot(rows,cols,2*i+1)
        nx.draw_spring(G, node_color = 'red')

    for i,G in enumerate(recomb_graphs_to_plot):

        G = nx.Graph(G) # copy the graph so we don't change anything in results
        remove_nodes_with_degree_0(G)

        plt.subplot(rows,cols,2*(i+1))
        nx.draw_spring(G, node_color = 'blue')

    # graphs = rows*cols
    # for i in range(graphs):
    #     plt.subplot(rows,cols,i)
    #     #print i
    #     results = simulate_single_graph_and_output_incompatibility_graph(n =3, L = 100, theta = 10.0, draw_largest_connected_component = True)
    #     #print
    #     determine_S_of_connected_components(results,print_results=False,with_header=True)
    plt.draw()
    # sample  incompatible datasets from each experiment, and plot them

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
    return float(max(G.degree().values()))/len(G.edges())

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
