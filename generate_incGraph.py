#!/usr/bin/python2.7
'''Generate incompatability graphs from MS simulations.

Written to generate test-data for my students in the Advanced Bioinformatics
course  2018.

(C) Mathias C. Cronjager 01.03.2018'''


from incompatability_graph_experiments import ms_sim_theta, incompatability_graph
from networkx import adjacency_matrix
import numpy as np

'''Define model parameters'''
n = 50 # number of sequences
theta = 100.0
rho = 30.0
sites_max = 100000

def mat2str(S):
    try:
        rowStrings = [ '[%s]'%', '.join(map(str,row)) for row in S ]
        outStr = '[%s]'%',\n '.join(rowStrings)
        return outStr
    except TypeError:
        return '[]'

def lengths2breakpoints(lengths):
    positions = []
    accumulator = 0.0
    for l in lengths[:-1]:
        accumulator += l
        positions.append(accumulator)

    return positions

sim = ms_sim_theta(n, theta, N = 1, rho = rho, seed = 1729)
command = sim['metadata']['command']
#S = sim['experiments'][0][0]
#loci = sim['experiments'][0][1]
lengths = [x['sites']/float(sites_max) for x in sim['experiments'][0][2]]
breakpoints = lengths2breakpoints(lengths)
blocks = len(lengths)

G = incompatability_graph(sim['experiments'][0])
loci   = sorted(G.nodes())

print 'command: %s\n'%command
print 'Sequences: %i'%n
print 'Mutation_rate: theta/2 = %g'%theta
print 'Recombination rate: rho/2 = %g\n'%rho
print 'Segregating sites: %i'%len(loci)
print 'No. tree-changes: %i'%(len(lengths) - 1)
print 'No. incompatible edge-pairs: %i'%len(G.edges())
print 'Positions (mutations):'
print loci
#print ', '.join(map(str,loci))
print '\nPositions (tree-changes):'
#print ', '.join(map(str,breakpoints))
print breakpoints
print '\nIncompatible position-pairs:'
print sorted(map(tuple,map(sorted,G.edges())))

# for row in S:
#     print row
