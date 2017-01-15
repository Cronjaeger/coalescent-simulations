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
    vertices = range(simulation.L)

    edges = simulation.getInconsistentColumnPairs()
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
    Will generate and output the characrters of all columns that feature in a connected components
    of an experiment. Each component is treated separately.

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
    p = subprocess.Popen(args,Shell = False, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    return output

def random_seed():
    '''
    an auxiliary function which allows us to generate a random seed to opass to
    ms, but which depends on the state of the ransom number generator used by
    np.random. This way, seeding the random number generator at the beginning of
    a script with "np.random.set_state()" will also guarantee that the seeds
    used in any sub-processes executring ms will be the same.
    '''
    return int(np.random.sample(1)*(2**31))

def ms_sim_theta(n, theta, N = 1, seed = random_seed(), extra_args = []):
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-t',str(theta)]+extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)

def ms_sim_sites(n, s, N = 1, seed = random_seed(), extra_args = []):
    args = ['ms',str(n),str(N),'-seeds',str(seed),'-s',str(s)] + extra_args
    output_raw = run_sim(args)
    output_parsed = parse_ms_output(output_raw)

def parse_ms_output(raw_input):
    'TODO: refine this!'
    return split(str(raw_input),'\n')

# def ms_sim(arg_string):
#     command = 'ms' + arg_string
#     p = subprocess.Popen(command,Shell = True, stdout=subprocess.PIPE)
#     output = p.communicate()[0]
#     #p.kill()
#     return output
