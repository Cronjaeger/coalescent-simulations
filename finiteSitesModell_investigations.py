# -*- coding: utf-8 -*-
"""
Created on Wed May 27 22:17:25 2015

@author: mathias
"""

import libCoal
import numpy as np
import matplotlib.pyplot as pl
from scipy.special import binom
from copy import deepcopy
from time import ctime
from math import log10
import networkx as nx

class coalescent_finiteSites(libCoal.coalescent):
    '''
    functionally identical to libCoal.coalescent. Only practical difference
    being that the site frequency spectrum is not computed.
    '''
    def computeSFS(self):
        pass

class simulator_KingmanFiniteSites(libCoal.simulateKingman):

    def __init__(self,n,mutationRate,L,compute_ancestral_configuration_on_initialization = True,*args):
        '''
        n = number of individuals
        mutationRate = Mutation rate of the coalescent.
        T_Max = Time-horizon of the coalescent.
        '''
        self.n = n
        self.mutationRate = mutationRate
#        self.mergerRate = mergerRate
        self.T_max = float('inf')
        self.T_MRCA = float('inf')
        self.args = args
        self.SFS = np.zeros(n)
        self.coal = coalescent_finiteSites(libCoal.partition([[i] for i in range(n)]))
#        self.simulateCoalescent(self.args)
        self.L = L
        self.MRCA_seq = np.zeros(self.L,dtype=int)
        self.site_mutationCount = np.zeros(self.L,dtype=int)
        #self.mutations_with_site_and_type = [] # lists of tuples of the form [(time,lineage,site,k)]
        self.sequences_mutationCount = np.zeros(self.n, dtype=int)
        self.sequences = np.zeros((self.n,self.L),dtype=int)
        self.compute_ancestral_configuration_on_initialization = compute_ancestral_configuration_on_initialization
        self.simulateCoalescent()

#    def preSimulationSteps(self):
#        self.L = self.args[0]
#        self.MRCA_seq = np.zeros(self.L,dtype=int)
#        self.site_mutationCount = np.zeros(self.L,dtype=int)
#        self.sequences_mutationCount = np.zeros(self.n, dtype=int)
#        self.sequences = np.zeros((self.n,self.L),dtype=int)

    # def getMutations_withSitesAndShifts(self):
    #     return deepcopy(self.mutations_with_site_and_type)

    def postSimulationSteps(self):
        if self.compute_ancestral_configuration_on_initialization:
            self.computeAncestralConfiguration()

    def computeAncestralConfiguration(self):
        affectedSites = np.random.randint(self.L, size=len(self.coal.mutations))
        mutationType = np.random.randint(1,4,size = len(self.coal.mutations))
        self.coal.mutations = [self.coal.mutations[i] + (affectedSites[i], mutationType[i],i) for i in xrange(len(self.coal.mutations))]
        for i in range(len(self.coal.mutations)):
            j = affectedSites[i]
            k = mutationType[i]
            self.site_mutationCount[j] += 1
            t, branch = self.coal.mutations[i][:2]
            affectedSequences = self.coal.getStateNoCoppy(t).blocks[branch]
            for sequenceIndex in affectedSequences:
                self.sequences[sequenceIndex,j] += k
                self.sequences[sequenceIndex,j] %= 4

                self.sequences_mutationCount[sequenceIndex] += 1

    def untillFirstXInconsistencies(self,computeS_seq = False,X = 2):

        """
        Step 1, set up everything like we were running
        computeAncestralConfiguration
        """
        M = len(self.coal.mutations)
        S = np.zeros((self.n,self.L),dtype=int)
        S_old = np.array(S)
        site_mut_count = np.zeros(self.L,dtype=int)
        seq_mut_count = np.zeros(self.n, dtype=int)
#        affectedSites = np.random.randint(self.L, size = M)
#        mutationType = np.random.randint(1,4,size = M)

        '''
        Step 2, randomly iterate over the list of mutations.
        Mutaions, untill the number of added mutations exceeds the number of
        segregating sites by X.
        '''
        inconsistencyCount = 0
        mutationCounter = 0
        ignoredMutations = list(self.coal.mutations)
        consideredMutations = []
        deviant_mutations = []
#        deviant_mutations_indices = []
#        self.coal.mutations = []
        typeCount =  [0,0,0,0]
        typeCountList = [list(typeCount)]
        S_seq = [np.matrix(S)]
#        inconsistentPairs = []
        """
        type 0 : a column with 3 states
        type 1 : a column with 2 states and >2 mutations
               : (no incompatibility w. other states)
        type 2 : creating a column with 2 states and incompatibilities
        type 3 : create an invisible state
        """
        while inconsistencyCount < X and mutationCounter < M:

            m_index = np.random.randint(M - mutationCounter)
            m_k = np.random.randint(1,4)
            m_site = np.random.randint(self.L)

            site_mut_count[m_site] += 1

            mutation = ignoredMutations.pop(m_index) + (m_site, m_k,mutationCounter)
            """
            Now a mutatiin has the following as its entries:
            m[0] : time of mutation
            m[1] : affected leneage
            m[2] : affected site
            m[3] : mutation type (k - shift mod 4)
            m[4] : a counter -- number of mutations preceeding this one
            """
            t,branch = mutation[:2]
            affectedSequences = self.coal.getStateNoCoppy(t).blocks[branch]

            for sequenceIndex in affectedSequences:
                S[sequenceIndex,m_site] += m_k
                S[sequenceIndex,m_site] %= 4
                seq_mut_count[sequenceIndex] += 1

#            if site_mut_count[m_site] > 1 and S[sequenceIndex,m_site] != 0:
            if site_mut_count[m_site] > 1:

                inconsistencyCount += 1
                mutation += (inconsistencyCount,)

                # inconsistent_columns_new = inconsistentColumnPairs(site_mut_count,S)
                # inconsistent_columns_old = inconsistentColumnPairs([ i - int(i==m_site) for i in site_mut_count],S_old)
                inconsistent_columns_new = inconsistentColumnPairs(S)
                inconsistent_columns_old = inconsistentColumnPairs(S_old)


                newInconsistencies = len(inconsistent_columns_new) > len(inconsistent_columns_old)

                if len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) > 1:
                    typeCount[0] += 1

                elif len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) == 1:
                    if newInconsistencies:
                        typeCount[2] += 1
                    else:
                        typeCount[1] += 1

                #if S[sequenceIndex,m_site] == 0:
                elif len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) == 0:
                    typeCount[3] += 1
                else:
                    raise Error("An exception which defies categorization has occurred. This should not happen!")

                typeCountList.append(list(typeCount))

                deviant_mutations.append(mutation)
#                deviant_mutations_indices.append(mutationCounter)

            mutationCounter += 1
            consideredMutations.append(mutation)
            if computeS_seq: S_seq.append(np.matrix(S))
            S_old[:,m_site] = S[:,m_site]

        if computeS_seq:
            newSeq = zip(consideredMutations,S_seq[1:])
            newSeq.sort(cmp = lambda x,y: int(np.sign(x[0][0] - y[0][0])))
            consideredMutations = [x[0] for x in newSeq]
            S_seq = S_seq[:1] + [x[1] for x in newSeq]
        else:
            consideredMutations.sort(cmp = lambda x,y: int(np.sign(x[0] - y[0])))
            newSeq = []

        self.coal.mutations = consideredMutations

        return {"S":S,
                "mutCount_sites":site_mut_count,
                "mutCount_sequences":seq_mut_count,
                "Inconsistencies":inconsistencyCount,
                "typeCount":typeCount,
                "typeCount_arr":np.array(typeCount),
                "coalescent" : deepcopy(self.coal),
                "typeCountList" : typeCountList,
                "newSeq" : newSeq,
                "S_seq" : S_seq,
                "deviant_mutations": deviant_mutations}

    def until_k_visible_mutations(self,k = 10):

        """
        Step 1, set up everything like we were running
        computeAncestralConfiguration
        """
        computeS_seq = False,
        M = len(self.coal.mutations)
        S = np.zeros((self.n,self.L),dtype=int)
        S_old = np.array(S)
        site_mut_count = np.zeros(self.L,dtype=int)
        seq_mut_count = np.zeros(self.n, dtype=int)
#        affectedSites = np.random.randint(self.L, size = M)
#        mutationType = np.random.randint(1,4,size = M)

        '''
        Step 2, randomly iterate over the list of mutations,
        untill enough mutations have occurred.
        '''
        inconsistencyCount = 0
        mutationCounter = 0
        ignoredMutations = list(self.coal.mutations)
        consideredMutations = []
        deviant_mutations = []
#        deviant_mutations_indices = []
#        self.coal.mutations = []
        typeCount =  [0,0,0,0]
        # typeCountList = [list(typeCount)]
        # S_seq = [np.matrix(S)]
#        inconsistentPairs = []
        """
        type 0 : a column with 3 states
        type 1 : a column with 2 states and >2 mutations
               : (no incompatibility w. other states)
        type 2 : creating a column with 2 states and incompatibilities
        type 3 : create an invisible state
        """
        while mutationCounter - (typeCount[1]+typeCount[3]) < k and mutationCounter < M:

            m_index = np.random.randint(M - mutationCounter)
            m_k = np.random.randint(1,4)
            m_site = np.random.randint(self.L)

            site_mut_count[m_site] += 1

            mutation = ignoredMutations.pop(m_index) + (m_site, m_k,mutationCounter)
            """
            Now a mutatiin has the following as its entries:
            m[0] : time of mutation
            m[1] : affected leneage
            m[2] : affected site
            m[3] : mutation type (k - shift mod 4)
            m[4] : a counter -- number of mutations preceeding this one
            """
            t,branch = mutation[:2]
            affectedSequences = self.coal.getStateNoCoppy(t).blocks[branch]

            for sequenceIndex in affectedSequences:
                S[sequenceIndex,m_site] += m_k
                S[sequenceIndex,m_site] %= 4
                seq_mut_count[sequenceIndex] += 1

#            if site_mut_count[m_site] > 1 and S[sequenceIndex,m_site] != 0:
            if site_mut_count[m_site] > 1:

                inconsistencyCount += 1
                mutation += (inconsistencyCount,)

                # inconsistent_columns_new = inconsistentColumnPairs(site_mut_count,S)
                # inconsistent_columns_old = inconsistentColumnPairs([ i - int(i==m_site) for i in site_mut_count],S_old)
                inconsistent_columns_new = inconsistentColumnPairs(S)
                inconsistent_columns_old = inconsistentColumnPairs(S_old)


                newInconsistencies = len(inconsistent_columns_new) > len(inconsistent_columns_old)

                if len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) > 1:
                    typeCount[0] += 1

                elif len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) == 1:
                    if newInconsistencies:
                        typeCount[2] += 1
                    else:
                        typeCount[1] += 1

                elif len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) == 0:
                    typeCount[3] += 1
                else:
                    raise Error("An exception which defies categorization has occurred. This should not happen!")

                # typeCountList.append(list(typeCount))

                deviant_mutations.append(mutation)
#                deviant_mutations_indices.append(mutationCounter)

            mutationCounter += 1
            consideredMutations.append(mutation)
            # if computeS_seq: S_seq.append(np.matrix(S))
            S_old[:,m_site] = S[:,m_site]

        # if computeS_seq:
        #     newSeq = zip(consideredMutations,S_seq[1:])
        #     newSeq.sort(cmp = lambda x,y: int(np.sign(x[0][0] - y[0][0])))
        #     consideredMutations = [x[0] for x in newSeq]
        #     S_seq = S_seq[:1] + [x[1] for x in newSeq]
        # else:
        #    consideredMutations.sort(cmp = lambda x,y: int(np.sign(x[0] - y[0])))
        #    newSeq = []
            consideredMutations.sort(cmp = lambda x,y: int(np.sign(x[0] - y[0])))
            # newSeq = []


        self.coal.mutations = consideredMutations

        return {"S":S,
                "mutCount_sites":site_mut_count,
                "mutCount_sequences":seq_mut_count,
                "Inconsistencies":inconsistencyCount,
                "typeCount":typeCount,
                "typeCount_arr":np.array(typeCount),
                "coalescent" : deepcopy(self.coal),
                # "typeCountList" : typeCountList,
                # "newSeq" : newSeq,
                # "S_seq" : S_seq,
                "deviant_mutations": deviant_mutations}

    def until_k_mutations(self,k = 10):
        '''
        Samples K mutations from self.muitations, and discards the rest.
        It then proceeds to compute the S-matrix associated with the coalescent.
        '''


        """
        Step 1, set up everything like we were running
        computeAncestralConfiguration
        """
        M = len(self.coal.mutations)
        S = np.zeros((self.n,self.L),dtype=int)
        S_old = np.array(S)
        # site_mut_count = np.zeros(self.L,dtype=int)
        # seq_mut_count = np.zeros(self.n, dtype=int)
#        affectedSites = np.random.randint(self.L, size = M)
#        mutationType = np.random.randint(1,4,size = M)
#
        '''
        Step 2, randomly iterate over the list of mutations.
        Mutaions, untill the number of added mutations equals k
        '''
        # inconsistencyCount = 0
        mutationCounter = 0
        ignoredMutations = list(self.coal.mutations)
        consideredMutations = []
        # deviant_mutations = []
#        deviant_mutations_indices = []
#        self.coal.mutations = []
        # typeCount =  [0,0,0,0]
        # typeCountList = [list(typeCount)]
        S_seq = [np.matrix(S)]
#        inconsistentPairs = []
        """
        type 0 : a column with 3 states
        type 1 : a column with 2 states and >2 mutations
               : (no incompatibility w. other states)
        type 2 : creating a column with 2 states and incompatibilities
        type 3 : create an invisible state
        """
        while mutationCounter <  min(M,k):

            m_index = np.random.randint(M - mutationCounter)
            m_k = np.random.randint(1,4)
            m_site = np.random.randint(self.L)

            # site_mut_count[m_site] += 1

            mutation = ignoredMutations.pop(m_index) + (m_site, m_k,mutationCounter)
            """
            Now a mutatiin has the following as its entries:
            m[0] : time of mutation
            m[1] : affected leneage
            m[2] : affected site
            m[3] : mutation type (k - shift mod 4)
            m[4] : a counter -- number of mutations preceeding this one
            """
            t,branch = mutation[:2]
            affectedSequences = self.coal.getStateNoCoppy(t).blocks[branch]

            for sequenceIndex in affectedSequences:
                S[sequenceIndex,m_site] += m_k
                S[sequenceIndex,m_site] %= 4
                # seq_mut_count[sequenceIndex] += 1

#            if site_mut_count[m_site] > 1 and S[sequenceIndex,m_site] != 0:
#             if site_mut_count[m_site] > 1:
#
#                 inconsistencyCount += 1
#                 mutation += (inconsistencyCount,)
#
#                 inconsistent_columns_new = inconsistentColumnPairs(site_mut_count,S)
#                 inconsistent_columns_old = inconsistentColumnPairs([ i - int(i==m_site) for i in site_mut_count],S_old)
#
#                 newInconsistencies = len(inconsistent_columns_new) > len(inconsistent_columns_old)
#
#                 if len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) > 1:
#                     typeCount[0] += 1
#
#                 if len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) == 1:
#                     if newInconsistencies:
#                         typeCount[2] += 1
#                     else:
#                         typeCount[1] += 1
#
#                 if S[sequenceIndex,m_site] == 0:
#                     typeCount[3] += 1
#
#                 typeCountList.append(list(typeCount))
#
#                 deviant_mutations.append(mutation)
# #                deviant_mutations_indices.append(mutationCounter)

            mutationCounter += 1
            consideredMutations.append(mutation)
            # if computeS_seq: S_seq.append(np.matrix(S))
            S_old[:,m_site] = S[:,m_site]

        # if computeS_seq:
        #     newSeq = zip(consideredMutations,S_seq[1:])
        #     newSeq.sort(cmp = lambda x,y: int(np.sign(x[0][0] - y[0][0])))
        #     consideredMutations = [x[0] for x in newSeq]
        #     S_seq = S_seq[:1] + [x[1] for x in newSeq]
        # else:
            # consideredMutations.sort(cmp = lambda x,y: int(np.sign(x[0] - y[0])))
            # newSeq = []
        consideredMutations.sort(cmp = lambda x,y: int(np.sign(x[0] - y[0])))
        #newSeq = []


        self.coal.mutations = consideredMutations
        self.sequences = S
#
        return {"S":S,
        #                 "mutCount_sites":site_mut_count,
        #                 "mutCount_sequences":seq_mut_count,
        #                 "Inconsistencies":inconsistencyCount,
        #                 "typeCount":typeCount,
        #                 "typeCount_arr":np.array(typeCount),
             "coalescent" : deepcopy(self.coal)
        #                 "typeCountList" : typeCountList,
        #                 "newSeq" : newSeq,
        #                 "S_seq" : S_seq,
        #                 "deviant_mutations": deviant_mutations
            }

    def getS(self):
        return np.array(self.sequences, dtype=int)

    def getSiteMutationCounts(self):
        return np.array(self.site_mutationCount, dtype = int)

    def getSeq_mutationCount(self):
        return np.array(self.sequences_mutationCount, dtype = int )

    def getSegregatingSites(self):
        return [j for j in range(self.L) if max(self.sequences[:,j]) > 0]

    def getInvisibleSites(self):
        return [j for j in range(self.L) if max(self.sequences[:,j])==0 and self.site_mutationCount[j] > 0]

    def countSegregatingSites(self):
#        counter = np.zeros(self.n, dtype = int)
        counter = 0
        for j in range(self.L):
            if len(filter(lambda x: x!= 0 , self.sequences[:,j])) > 0:
                counter += 1
        return counter

    def countMinimalMutations(self):
        counter = 0
        for j in range(self.L):
            mutations = filter(lambda x: x!= 0 , self.sequences[:,j])
            if len(mutations) > 0:
                counter += len(set(mutations))
        return counter

    def countInconsistencies(self):
        #return len(inconsistentColumnPairs(site_mut_count = self.site_mutationCount , S = self.sequences) )
        return len(inconsistentColumnPairs(S = self.sequences) )

    def getInconsistentColumnPairs(self):
        #return inconsistentColumnPairs(site_mut_count = self.site_mutationCount , S = self.sequences)
        return inconsistentColumnPairs(S = self.sequences)

# def inconsistentColumnPairs(site_mut_count, S, ancestral_type_known = True):
#     if S.shape[1] != len(site_mut_count):
#         raise ValueError('Number of collumns in S (=%i) does not match length of site_mut_count (=%i)'%(S.shape[1],len(site_mut_count)))
#     pairs = []
#     affectedSites = filter(lambda i: sum(S[j,i] != 0 for j in xrange(S.shape[0])) > 1, xrange(S.shape[1]) )
#     #affectedSites = filter(lambda i: site_mut_count[i] > 0, xrange(len(site_mut_count)) )
#     for s1 in affectedSites:
#         for s2 in filter(lambda x: x > s1 , affectedSites):
#             #if not three_gammete_test(S[:,s1],S[:,s2]):
#             if not compatibility_test(S[:,s1], S[:,s2], ancestral_type_known = ancestral_type_known):
#                 pairs.append((s1,s2))
#     return pairs

def inconsistentColumnPairs(S, ancestral_type_known = True):
    pairs = []
    number_of_non_0_elements = np.apply_along_axis(func1d = sum, axis = 0, arr = (S != 0))
    affectedSites = [i for i,n in enumerate(number_of_non_0_elements) if n > 1]
    # discriminant = lambda i: sum(S[j,i] != 0 for j in xrange(S.shape[0])) > 1
    # affectedSites = filter(discriminant, xrange(S.shape[1]))
    for i,s1 in enumerate(affectedSites):
        #for s2 in filter(lambda x: x > s1 , affectedSites):
        for s2 in affectedSites[(i+1):]:
            if not two_char_compatibility_test(S[:,s1], S[:,s2], ancestral_type_known = ancestral_type_known):
                pairs.append((s1,s2))
    return pairs


def two_char_compatibility_test(c1,c2, ancestral_type_known = False, ancestral_type = (0,0)):
    '''
    Takes two columns of equal length (representing characters on a phylogeny),
    and returns true if and only if they are compatible (in the sense that the
    number of mutations nessecary to place both characters on a single tree is
    the sum of the numbr of mutations nessecary to place each character on a
    tree individually)
    '''
    if len(c1) != len(c2): raise(ValueError('Cannot compare characters of unequal length'))

    #convert input to list
    c1_list = list(c1)
    c2_list = list(c2)

    if ancestral_type_known:
        c1_list.append(ancestral_type[0])
        c2_list.append(ancestral_type[1])

    # G_dict = partition_intersection_graph([c1_list,c2_list])
    # G_nx = nx.Graph(G_dict['E'])
    G_nx = partition_intersection_graph([c1_list,c2_list])

    return nx.is_forest(G_nx)
    #return len(nx.cycles.cycle_basis(G_nx)) == 0


def three_gammete_test(c1,c2):
    '''
    Takes two columns/or characters.
    Returns False if the characters contain all three of the gammetes (not 0, 0)
    (not 0, not 0) and (0, not 0). Returns True otherwise.
    '''
    n = len(c1)

    AX_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] == 0 and c2[i] != 0, range(n)))
    YA_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] != 0 and c2[i] == 0, range(n)))
    YX_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] != 0 and c2[i] != 0, range(n)))

    return not(AX_occurs and YA_occurs and YX_occurs)

def compatibility_criterion_1(c1,c2, ancestral_type_known = False, ancestral_type = (0,0)):
    '''
    THIS CRITERION IS NESSECARY BUT INSUFFICIENT FOR PAIRWISE COMPATIBILITY!
    Takes two columns of equal length (representing characters on a phylogeny),
    and returns true if (NOT only if) they are compatible (in the sense that the
    number of mutations nessecary to place both characters on a single tree is
    the sum of the numbr of mutations nessecary to place each character on a
    tree individually)
    '''
    if len(c1) != len(c2): raise(ValueError('Cannot compare characters of unequal length'))

    states_1 = set(c1)
    states_2 = set(c2)

    gammetes = set(zip(c1,c2))

    if ancestral_type_known:
        gammetes.add((0,0)) # gammetes is a set; only has effect when (0,0) not already present.

    n_states_1 = len(states_1)
    n_states_2 = len(states_2)

    n_gammetes = len(set(gammetes))

    return n_gammetes <= n_states_1 + n_states_2 - 1 #,states_1,states_2,gammetes

def partition_intersection_graph(chars, output_as_dict = False):
    '''
    accepts a collection (e.g. a list) of characters, and generates their
    partition intersection graph G. The output is either a dictionary, where
    G['V'] = nodes, and G['E'] = edges, or an object of the type networkx.Graph.
    '''
    nodes = []
    for i,char in enumerate(chars):
        p = to_partition(char)
        for block in p:
            nodes.append((i,tuple(block)))

    edges = []
    for i, a in enumerate(nodes):
        for b in nodes[i+1:]:
            if len(set(a[1]).intersection(set(b[1]))) > 0:
                edges.append([a,b])

    if output_as_dict:
        return {'V':nodes,'E':edges}
    else:
        return nx.Graph(edges)

def get_informative_collumns(S):
    '''A columns where all but one character state exists only once is deemed
    'un-informative'. This method removes all such columns.
    '''
    number_of_non_0_elements = np.apply_along_axis(func1d = sum, axis = 0, arr = (S != 0))
    segregatingSites = [i for i,n in enumerate(number_of_non_0_elements) if n>1]
    #return filter(lambda i: is_informative(S[:,i]), range(S.shape[1]))
    return filter(lambda i: is_informative(S[:,i]), segregatingSites)

def get_non_informative_collumns(S):
    '''
    returns the complement of get_informative_collumns(S)
    '''
    return filter(lambda i: not is_informative(S[:,i]), range(S.shape[1]))

def get_essential_state_counts(S):
    return map(lambda i: number_of_essential_states(S[:,i]), range(S.shape[1]))

def is_informative(char,ancestral_type_known = True, ancestral_type = 0):
    ''' Determine if a character is informative

    Returns: Boolean value
        True iff char has only one block, or multiple blocks, with the
        seccond-largest block having size 1 (i.e. if only one non-trivial block)
    '''
    char = list(char)
    if ancestral_type_known:
        char.append(ancestral_type)

    n_part = map(len,to_partition(char)) # get block-sizes of partition
    if len(n_part) == 1:
        return False # len(part)==1 iff single state character
    else:
        n_part.sort(reverse = True)
        return n_part[1] > 1

def number_of_essential_states(char, ancestral_type_known = True, ancestral_type = 0):
    '''returns the number of states which appear more than one.
if the ancestral type is indicated as known, the count of the state indicated by
the variable 'ancestral_type' artificially has its count boosted by 1.
    '''
    char = list(char)
    if ancestral_type_known:
        char = char + [ancestral_type]
    n_part = map(len,to_partition(char))
    return sum([b > 1 for b in n_part])


def to_partition(char, block_constructor = set):
    '''in : a character e.g. ['A','A','G','A','T']
out: a list of sets. each set corresponds to all indices (from 0 to len(char)-1)
with the same state, e.g. [set([0,1,3]),set([2]),set([4])]. This list will be
sorted in order of least elements.'''
    states = set(list(char))
    partition = [block_constructor( filter( lambda i: char[i] == s, range(len(char)) ) ) for s in states]
    partition.sort(key = min)
    return partition

def count_minimal_vs_actual_mutations(simulation, ancestral_type_known = True, ret_all= True):
    #assert isinstance(simulation,simulator_KingmanFiniteSites)

    site_mut_count = simulation.getSiteMutationCounts() # counts actual number of mutations

    S = simulation.getS()
    if ancestral_type_known:
        # add a (0,0,0,...,0)-row to S; representing the sequence of the MRCA
        S = np.r_[S,np.zeros( (1,S.shape[1]), dtype = int)]

    block_counts = [len( to_partition( S[:,i] ) ) for i in range(S.shape[1])]
    min_mut_counts = [i-1 for i in block_counts]

    A = np.array(site_mut_count)
    B = np.array(min_mut_counts)
    C = A - B
    if ret_all:
        return {'actual':A, 'minimal':B, 'diff':C}
    else:
        return C


def draw_partition_intersection_graph(representation):
    if type(representation) == list: # we have been passed a list of characters
        chars = representation
        G = partition_intersection_graph(chars)
    else:
        assert type(representation) == nx.Graph
        G = representation

    nx.draw(G, node_color = [n[0] for n in G.nodes()])

def fromEdgesToConnectedComponents(pairs):
    '''
    In: a list of edges (i,j) satisfying i < j
    Out: a list of sets [S_1, S_2, ...] corresponding to the verticees of
connected components in the graph encoded by the input.
    '''
    if len(pairs) == 0:
        return []

    # vertices_with_duplicates = reduce( lambda x,y: x+y, pairs)
    # verticees = list(set(vertices_with_duplicates))
    # verticees.sort()


    #we encode the set of all blocks in the graph by the smallest vertex if.
    # hence j belongs to the block encoded by i if and only if is the minimum
    # index in the connected component containing j.
    blocks = {}
    block_index = {}
    # blocks = dict(verticees[0]:set(verticees[0],))
    # block_indices = dict(verticees[0]:verticees[0])
    #for v in verticees:
    # verticees_seen = set([])
    for edge in pairs:

        vertices_seen = block_index.keys()

        v1 = edge[0]
        v2 = edge[1]

        if v1 not in vertices_seen and v2 not in vertices_seen:
            block_index[v1] = v1
            block_index[v2] = v1
            blocks[v1] = set([v1,v2])

        if v1 in vertices_seen and v2 not in vertices_seen:
            block_index[v2] = block_index[v1]
            blocks[block_index[v1]].add(v2)

        if v1 not in vertices_seen and v2 in vertices_seen:
            block_index[v1] = block_index[v2]
            blocks[block_index[v2]].add(v1)

        if v1 in vertices_seen and v2 in vertices_seen:
            if block_index[v1] != block_index[v2]:

                if block_index[v1] < block_index[v2]:
                    index_min = block_index[v1]
                    index_max = block_index[v2]
                else:
                    index_max = block_index[v1]
                    index_min = block_index[v2]

                block_max = blocks.pop(index_max)
                blocks[index_min] = blocks[index_min].union(block_max)

                for v in block_max:
                    block_index[v] = index_min

    return blocks.values()

def printPairListWithDegrees(pairs,components):
    vertices_with_duplicates = reduce( lambda x,y: x+y, pairs)
    verticees = list(set(vertices_with_duplicates))
    verticees.sort()

    degree = dict([(v,0) for v in verticees])

    for v1,v2 in pairs:
        degree[v1] += 1
        degree[v2] += 1

    degreeList = degree.items()
    degreeList.sort(cmp = lambda x,y: y[1] - x[1])

    print 'vertices (%i total):'%len(verticees)
    print verticees
    print '\nedges (%i total):'%len(pairs)
    print set(pairs)
    print '\nconnected components:(%i total)'%len(components)
    print components
    print '\ndegrees (as (vertex,degree) ):'
    print degreeList

def generate_plots_for_jotun(n = 8,L = 100,thetaMax = 100,thetaMin=0.001,steps=5,N=100,savePath = ''):

    saveFigures = len(savePath) > 0

    print_a_larger_than_10 = 10
    print_a_50 = True

    #Run simulations
    thetas = np.logspace(log10(thetaMin),log10(thetaMax),steps)
    typeCouunter_simulationAverages = np.zeros((len(thetas),4))
    averageInconsistencyBlocksizeFrequency = np.zeros((len(thetas),L))
    for i,theta in enumerate(thetas):

        simulations = [simulator_KingmanFiniteSites(n,float(theta)/2,L) for z in range(N)]
        typeCounter = np.zeros((N,4))

        for j in range(N):

            s = simulations[j]
            S = s.getS()

            for column in range(S.shape[1]):
                typeCounter[j,len(set(S[:,column]))-1] += 1


            incompatible_pairs = s.getInconsistentColumnPairs()
            incompatible_components = fromEdgesToConnectedComponents(incompatible_pairs)
            # print incompatible_pairs
            # print incompatible_components
            incompatability_group_sizes = map(len,incompatible_components)
            #incompatability_group_sizes = [len(c) for c in incompatible_components]


            for k in incompatability_group_sizes:
                if print_a_larger_than_10 > 0 and k >= 10:
                    print 'Incompatitility graph: generated with (n,L,theta) = (%i,%i,%10.5f):\n'%(n,L,theta)
                    printPairListWithDegrees(incompatible_pairs,incompatible_components)
                    print_a_larger_than_10 -= 1
                    print '-'*80
                averageInconsistencyBlocksizeFrequency[i][k -1] += 1

        typeCouunter_simulationAverages[i] = [sum(typeCounter[:,k]) / float(N) for k in range(4)]

        averageInconsistencyBlocksizeFrequency[i] *= 1.0/N

    #generate plot of frequency of number of dirrerent characters:
    fig = pl.figure(figsize=(20,5))
    pl.suptitle('%i sequences of length %i; %i simulations per point '%(n,L,N))
    pl.subplot(121)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel("mean frequency of loci with k types")
    pl.xscale('log')
    for k in (1,2,3,4):
        pl.plot(thetas/L , typeCouunter_simulationAverages[:,k-1], label = '$k = %i$'%k)
    pl.legend(loc='best')

    #generate plot of frequency of incompatability-blocksizes
    pl.subplot(122)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel("mean frequency of maximal\nincompatible groups of size s")
    pl.xscale('log')
    blocksizelist = [2,3,4,5,10,25,50]
    for k in blocksizelist:
        pl.plot(thetas/L , averageInconsistencyBlocksizeFrequency[:,k-1],label = '$s = %i$'%k)
    pl.legend(loc='best')

    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    if saveFigures:
        pl.savefig(savePath+"jotunPlot__L_%i__N_%i__n_%i.pdf"%(L,N,n))
        pl.savefig(savePath+"jotunPlot__L_%i__N_%i__n_%i.png"%(L,N,n))
        pl.savefig(savePath+"jotunPlot__L_%i__N_%i__n_%i.eps"%(L,N,n))
        pl.savefig(savePath+"jotunPlot__L_%i__N_%i__n_%i.svg"%(L,N,n))

        csv_path = savePath+"jotunPlot__L_%i__N_%i__n_%i"%(L,N,n)
        #csv_out = open(csv_path,w)
        array1_out = np.c_[thetas/L , typeCouunter_simulationAverages]
        np.savetxt(csv_path+'_subfig_1.csv',array1_out,fmt='%10.10f',delimiter = ', ')
        array2_out = np.c_[thetas/L , averageInconsistencyBlocksizeFrequency[:,blocksizelist]]
        np.savetxt(csv_path+'_subfig_2.csv',array2_out,fmt='%10.10f',delimiter = ', ')
        # for vector in [thetas/L]+[[]]+[typeCouunter_simulationAverages[:,k-1] for k in (1,2,3,4)]+[[]]+[averageInconsistencyBlocksizeFrequency[:,k-1] for k in range(2,11)]:
        #     csv_out.write(', '.join(['%.10f'%x for x in vector]))
        #     csv_out.write('\n')
        # csv_out.close()
    pl.show()

def generate_plot_1(n = 8,L = 100,thetaMax = 10,thetaMin=0.01,steps=20,N=100,savePath = ''):

    saveFigures = len(savePath) > 0

    #Run simulations
    h = (float(thetaMax) - thetaMin)/steps
    thetas = np.arange(thetaMin,thetaMax,h)+h
    avgRate = np.zeros(len(thetas))
    inconsistencyCount = np.zeros(len(thetas))
    invisibleSiteCount = np.zeros(len(thetas))
    for i,theta in enumerate(thetas):
        simulations = [simulator_KingmanFiniteSites(n,float(theta)/2,L) for z in range(N)]
        # rates = []
        # inconsistencies = []
        # invisibleSites = 0
        # segCounter = 0
        for s in simulations:

#        for i in range(N):
#            s = simulations[i]

            minimalMutations = s.countMinimalMutations()
            actualMutations = len(s.coal.mutations)
            segregatingSites = s.countSegregatingSites()

#            siteMutationCounts = s.getSiteMutationCounts()

            if minimalMutations > 0:
                rates.append( float(actualMutations) / minimalMutations )

                inconsistencies.append(s.countInconsistencies())

                invisibleSites += len(s.getInvisibleSites())/float(segregatingSites)
                segCounter += 1

        if segCounter > 0:
            invisibleSiteCount[i] = invisibleSites/float(segCounter)
        else:
            invisibleSiteCount[i] = 0

        avgRate[i] = np.average(rates)
        inconsistencyCount[i] = np.average(inconsistencies) / binom(L,2)

    #generate plot 1
    pl.figure()
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"(actual # mutations) / (# visible mutations)")
    pl.plot(thetas/L , avgRate , color='blue' , label=label)
    pl.legend(loc='upper left')
    if saveFigures: pl.savefig(savePath+"plot1__L_%i__N_%i__n_%i.pdf"%(L,N,n))

    #generate plot 2
    pl.figure()
    #    label = "(x;y) = (theta/L ; fraction of inconsistent columns)\nL,N,n = %i,%i,%i"%(L,N,n)
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"#inconsistent column-pairs / $\binom{L}{2}$")
    pl.plot(thetas/L, inconsistencyCount, color = "red", label = label)
    pl.legend(loc='upper left')
    if saveFigures: pl.savefig(savePath+"plots/plot2__L_%i__N_%i__n_%i.pdf"%(L,N,n))

    #generate plot 3
    pl.figure()
    #    label = "(x;y) = (theta/L ; fraction of invisible sites)\nL,N,n = %i,%i,%i"%(L,N,n)
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"#invisible sites / #segregating sites")
    pl.plot(thetas/L, invisibleSiteCount, color = "green", label = label)
    pl.legend(loc='upper right')
    if saveFigures: pl.savefig(savePath+"plots/plot3__L_%i__N_%i__n_%i.pdf"%(L,N,n))
    pl.show()

#def generate_plot_2(n,L,theta,N=100):
#    simulations = [simulator_KingmanFiniteSites(n,float(theta)/2,L) for z in range(N)]
#    count_inconsistencies_total = 0
#    for i,s in enumerate(simulations):
#        count_inconsistencies = s.countInconsistColumnPairs()
#        count_inconsistencies_total += count_inconsistencies


def runTests():
#    generate_plot_1(n=10,L=50,thetaMax=50,steps=50,N=1000)
    generate_plot_1(n=20,L=200,thetaMax=200,steps=20,N=1000)

def simulateUntillXMutations(N = 1000, n = 20, L = 100, mutRate = 200,printFirst10 = False,X = 2):
    K_list = [simulator_KingmanFiniteSites(n,mutRate,L,False) for i in xrange(N)]
    totalTypeCount = np.array((0,0,0,0))
    misses = 0
    k_res_List = []

    for K in K_list:
        res = K.untillFirstXInconsistencies(X = X)
        #Guarentee that we actually got enough mutations.
        if res["Inconsistencies"] == X:
            totalTypeCount += res["typeCount_arr"]
            k_res_List.append([K,res])
        else:
            misses += 1

    if printFirst10:
        for K in K_list[:10]:
            print chronology(K)
            print ""

    N_eff = N - misses

    return totalTypeCount,N_eff,k_res_List

def chronology(K):
    '''
    Outputs a String of all events in K, sorted in chronological order
    '''
    events = list(K.coal.jumps)
    events.extend(K.coal.mutations)
    events.sort(cmp = lambda x,y: int(np.sign(x[0] - y[0])))
    return "\n".join(map(eventToString,events))

def eventToString(e):
    if str(type(e[1])) == "<class 'libCoal.partition'>":
        return "t = %1.4f\t"%e[0] + str(e[1])
    if isinstance(e[1], int):
        outStr = "t = %1.4f\t"%e[0]
        if len(e) > 5: outStr += " *** Inconsistency number %i *** "%e[5]
        outStr += "Mutation (%i in order added"%e[4]
        outStr += "; +%i-mod4 at site %i affecting block %i)"%(e[3],e[2],e[1])
        return outStr
    else:
        print "WTF!",e

def scatterplot_index_and_time_of_abnormal_mutations(N = 1000 , L = 100, n = 20):

    theta = 1.2 * L
    typeCounts,N_eff,k_res_list = simulateUntillXMutations(N = N, n = n , L = L, mutRate=theta, printFirst10= False)
    t1, I1 = np.zeros(N_eff,dtype = float), np.zeros(N_eff,dtype = int)
    t2, I2 = np.zeros(N_eff,dtype = float), np.zeros(N_eff,dtype = int)
    i = 0
    for K,res in k_res_list:
        mutations =  res["deviant_mutations"]
        t1[i], I1[i] = mutations[0][0],mutations[0][4]
        t2[i], I2[i] = mutations[1][0],mutations[1][4]
        i += 1

    pl.figure()

    pl.subplot(1,2,1)
    pl.scatter(t1,t2,marker="+",s=10)
    pl.axis('scaled')
    pl.xlim(xmin = 0.0)
    pl.ylim(ymin = 0.0)
    pl.title("(coalescent) time of mutations")
    pl.xlabel(r"$t_1$")
    pl.ylabel(r"$t_2$")

    pl.subplot(1,2,2)
    pl.scatter(I1,I2,marker="+",s=10)
    pl.xlim(xmin = 0.0)
    pl.ylim(ymin = 0.0)
    pl.axis('scaled')
    pl.xlabel(r"$i_1$")
    pl.ylabel(r"$i_2$")
    pl.title("Index of mutations\n(in order of simulation)")

    pl.suptitle("N = %i      L = %i      n = %i"%(N_eff,L,n))

    filename_str = "plots/scatter_mut1_vs_2/scatter_mut1_vs_2_N_%i_L_%i_n_%i"%(N_eff,L,n)
    try:
        pl.savefig(filename_str+".pdf")
        pl.savefig(filename_str+".png")
        pl.savefig(filename_str+".ps")
        pl.savefig(filename_str+".svg")
        pl.savefig(filename_str+".eps")
    except Exception:
        print "could not save in all formats"

    pl.draw()


def generatePlot_of_mutationTypes(N = 1000,L = 100, n = 20, printFirsrst10 = False,show = False, X=2):

    theta = 2.0 * L
    typeCounts,N_eff,K_list = simulateUntillXMutations(N = N, n = n , L = L, mutRate=theta, printFirst10= printFirsrst10, X = X)

    #run simulations
    #typeCounts,N_eff = simulateUntillXMutations(N = N, n = n , L = L, mutRate=theta, printFirst10= printFirsrst10,X = X)[:2]

    #plot simulation-results
    width = 0.8
    #color = ("cyan","orange","grey","magenta")
    color = ("white","white","white","white")
    hatch = ("//", ".","\\")
    left = np.arange(1,5) - width/2.0
    pl.figure()
    bars = pl.bar(left,typeCounts,width = width, color = color)
    for bar,pattern in zip(bars[1:],hatch):
        bar.set_hatch(pattern)
#    pl.xlabel("Type of incompatibility")
    pl.xticks(np.arange(1,5),("3 types","2 types\n2 mutations\n3 gammete test passes","2 types\n3 gammete test fails","1 type\n2 mutations"))
    pl.ylabel("frequency")
    pl.title("Result of %i simulations stopped after %i events\nsequences = %i    sequence-length = %i"%(N_eff,X,n,L))
#    pl.tight_layout()
    pl.draw()
    filename_str = "plots/bars_stoppedProcess/bar_typeFrequencies_X_%i_N_%i_L_%i_n_%i_no_color"%(X,N_eff,L,n)
    try:
        pl.savefig(filename_str+".pdf")
        pl.savefig(filename_str+".png")
        pl.savefig(filename_str+".ps")
        pl.savefig(filename_str+".svg")
        pl.savefig(filename_str+".eps")
    except Exception:
        print "could not save in all formats (pdf,png,ps,svg,eps)"
    if show:
        pl.show()
    np.savetxt(filename_str+".csv", np.array(typeCounts,ndmin=2),fmt='%d', delimiter = ', ',header = "(a) 3 types, (b) 2 types 2 mutations 3 gammete test passes, (c) 2 types 3 gammete test fails, (d) 1 type 2 mutations")


def run_generatePlot_of_mutationTypes(arglist = [(1000,100,20)],X = 2):
    for args in arglist:
        N,L,n = args
        print "%s   Generating plots (Barcharts) for N,L,n,X = %i,%i,%i,%i"%((ctime(),)+tuple(args)+(X,))
        generatePlot_of_mutationTypes(N=N,L=L,n=n,printFirsrst10 = False,show = False, X = X)

def run_generateScatterplots(arglist = [(100,100,20)]):
    for args in arglist:
        N,L,n = args
        print "Generating plots (scatterplts) for N,L,n = %i,%i,%i"%tuple(args)
        scatterplot_index_and_time_of_abnormal_mutations(N,L,n)
