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

class coalescent_no_SFS(libCoal.coalescent):
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
        self.coal = coalescent_no_SFS(libCoal.partition([[i] for i in range(n)]))
#        self.simulateCoalescent(self.args)
        self.L = L
        self.MRCA_seq = np.zeros(self.L,dtype=int)
        self.site_mutationCount = np.zeros(self.L,dtype=int)
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


    def postSimulationSteps(self):
        if self.compute_ancestral_configuration_on_initialization:
            self.computeAncestralConfiguration()

    def computeAncestralConfiguration(self):
        affectedSites = np.random.randint(self.L, size=len(self.coal.mutations))
        mutationType = np.random.randint(1,4,size = len(self.coal.mutations))
        for i in range(len(self.coal.mutations)):
            j = affectedSites[i]
            k = mutationType[i]
            self.site_mutationCount[j] += 1
            t,branch = self.coal.mutations[i]
            affectedSequences = self.coal.getStateNoCoppy(t).blocks[branch]
            for sequenceIndex in affectedSequences:
                self.sequences[sequenceIndex,j] += k
                self.sequences[sequenceIndex,j] %= 4

                self.sequences_mutationCount[sequenceIndex] += 1

    def untillFirstTwoInconsistencies(self):

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
        segregating sited by two.
        '''
        inconsistencyCount = 0
        mutationCounter = 0
        ignoredMutations = list(self.coal.mutations)
        typeCount =  [0,0,0,0]
#        inconsistentPairs = []
        """
        type 0 : a column with 3 states
        type 1 : a column with 2 states and >2 mutations
               : (no incompatibility w. other states)
        type 2 : creating a column with 2 states and incompatibilities
        type 3 : create an invisible state
        """
        while inconsistencyCount < 2 and mutationCounter < M:

            m_index = np.random.randint(M - mutationCounter)
            m_k = np.random.randint(1,4)
            m_site = np.random.randint(self.L)

            site_mut_count[m_site] += 1

            t,branch = ignoredMutations.pop(m_index)
            affectedSequences = self.coal.getStateNoCoppy(t).blocks[branch]

            for sequenceIndex in affectedSequences:
                S[sequenceIndex,m_site] += m_k
                S[sequenceIndex,m_site] %= 4
                seq_mut_count[sequenceIndex] += 1

#            if site_mut_count[m_site] > 1 and S[sequenceIndex,m_site] != 0:
            if site_mut_count[m_site] > 1:

                inconsistencyCount += 1

                inconsistent_columns_new = inconsistentColumnPairs(site_mut_count,S)
                inconsistent_columns_old = inconsistentColumnPairs([ i - int(i==m_site) for i in site_mut_count],S_old)

                newInconsistencies = len(inconsistent_columns_new) > len(inconsistent_columns_old)

                if len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) > 1:
                    typeCount[0] += 1

                if len(set((S[i,m_site] for i in xrange(S.shape[0]))) - set((0,))) == 1:
                    if newInconsistencies:
                        typeCount[2] += 1
                    else:
                        typeCount[1] += 1

                if S[sequenceIndex,m_site] == 0:
                    typeCount[3] += 1

            mutationCounter += 1
            S_old[:,m_site] = S[:,m_site]

        return {"S":S,
                "mutCount_sites":site_mut_count,
                "mutCount_sequences":seq_mut_count,
                "Inconsistencies":inconsistencyCount,
                "typeCount":typeCount,
                "typeCount_arr":np.array(typeCount)}

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
        return len(inconsistentColumnPairs(site_mut_count = self.site_mutationCount , S = self.sequences) )

def inconsistentColumnPairs(site_mut_count, S):
    pairs = []
    affectedSites = filter(lambda i: site_mut_count[i] > 0, xrange(len(site_mut_count)) )
    for s1 in affectedSites:
        for s2 in filter(lambda x: x > s1 , affectedSites):
            if isInconsistentColumnPair(S[:,s1],S[:,s2]):
                pairs.append((s1,s2))
    return pairs

def isInconsistentColumnPair(c1,c2):

    n = len(c1)

    AX_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] == 0 and c2[i] != 0, range(n)))
    YA_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] != 0 and c2[i] == 0, range(n)))
    YX_occurs = reduce(lambda x,y: x or y, map(lambda i: c1[i] != 0 and c2[i] != 0, range(n)))

    return (AX_occurs and YA_occurs and YX_occurs)


def generate_plot_1(n,L,thetaMax,thetaMin=0,steps=20,N=100):

    #Run simulations
    h = (float(thetaMax) - thetaMin)/steps
    thetas = np.arange(thetaMin,thetaMax,h)+h
    avgRate = np.zeros(len(thetas))
    inconsistencyCount = np.zeros(len(thetas))
    invisibleSiteCount = np.zeros(len(thetas))
    for i,theta in enumerate(thetas):
        simulations = [simulator_KingmanFiniteSites(n,float(theta)/2,L) for z in range(N)]
        rates = []
        inconsistencies = []
        invisibleSites = 0
        segCounter = 0
        for s in simulations:
#        for i in range(N):
#            s = simulations[i]

            minimalMutations = s.countMinimalMutations()
            actualMutations = len(s.coal.mutations)
            segregatingSites = s.countSegregatingSites()

#            siteMutationCounts = s.getSiteMutationCounts()

            if minimalMutations > 0:
                rates.append( float(actualMutations) / minimalMutations )

                inconsistencies.append(s.countInconsistColumnPairs())

                invisibleSites += len(s.getInvisibleSites())/float(segregatingSites)
                segCounter += 1

        invisibleSiteCount[i] = invisibleSites/float(segCounter)
        avgRate[i] = np.average(rates)
        inconsistencyCount[i] = np.average(inconsistencies) / binom(L,2)

    #generate plot 1
    pl.figure()
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"(actual # mutations) / (# visible mutations)")
    pl.plot(thetas/L , avgRate , color='blue' , label=label)
    pl.legend(loc='upper left')
    pl.savefig("plots/plot1__L_%i__N_%i__n_%i.pdf"%(L,N,n))

    #generate plot 2
    pl.figure()
#    label = "(x;y) = (theta/L ; fraction of inconsistent columns)\nL,N,n = %i,%i,%i"%(L,N,n)
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"#inconsistent column-pairs / $\binom{L}{2}$")
    pl.plot(thetas/L, inconsistencyCount, color = "red", label = label)
    pl.legend(loc='upper left')
    pl.savefig("plots/plot2__L_%i__N_%i__n_%i.pdf"%(L,N,n))

    #generate plot 3
    pl.figure()
#    label = "(x;y) = (theta/L ; fraction of invisible sites)\nL,N,n = %i,%i,%i"%(L,N,n)
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"#invisible sites / #segregating sites")
    pl.plot(thetas/L, invisibleSiteCount, color = "green", label = label)
    pl.legend(loc='upper right')
    pl.savefig("plots/plot3__L_%i__N_%i__n_%i.pdf"%(L,N,n))


#def generate_plot_2(n,L,theta,N=100):
#    simulations = [simulator_KingmanFiniteSites(n,float(theta)/2,L) for z in range(N)]
#    count_inconsistencies_total = 0
#    for i,s in enumerate(simulations):
#        count_inconsistencies = s.countInconsistColumnPairs()
#        count_inconsistencies_total += count_inconsistencies


def runTests():
#    generate_plot_1(n=10,L=50,thetaMax=50,steps=50,N=1000)
    generate_plot_1(n=20,L=200,thetaMax=200,steps=20,N=1000)

def simulateUntillTwoMutations(N = 1000, n = 20, L = 100, mutRate = 200):
    K_list = [simulator_KingmanFiniteSites(n,mutRate,L,False) for i in xrange(N)]
    totalTypeCount = np.array((0,0,0,0))
    misses = 0

    for K in K_list:
        res = K.untillFirstTwoInconsistencies()
        if res["Inconsistencies"] == 2:
            totalTypeCount += res["typeCount_arr"]
        else:
            misses += 1

    return totalTypeCount,N-misses

def generatePlot_of_mutationTypes(N = 1000,L = 100, n = 20):

    theta = 1.2 * L

    #run simulations
    typeCounts,N_eff = simulateUntillTwoMutations(N = N, n = n , L = L, mutRate=theta)

    #plot simulation-results
    width = 0.8
    color = ("cyan","orange","grey","magenta")
    left = np.arange(1,5) - width/2.0
    pl.figure()
    pl.bar(left,typeCounts,width = width, color = color)
#    pl.xlabel("Type of incompatibility")
    pl.xticks(np.arange(1,5),(">2 types","2 types\n>2 mutations\nno incompatibility","2 types\nincompatibility","invisible site"))
    pl.ylabel("frequency")
    pl.title("Result of %i simulations stopped after 2 events\nsequences = %i    sequence-length = %i"%(N_eff,n,L))
#    pl.tight_layout()
    pl.show()
    pl.savefig("plots/bars_stoppedProcess/bar_typeFrequencies_N_%i_L_%i_n_%i.pdf"%(N_eff,L,n))
