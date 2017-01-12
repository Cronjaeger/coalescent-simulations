import networkx as nx
import matplotlib.pyplot as plt
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
