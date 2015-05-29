# -*- coding: utf-8 -*-
"""
Created on Wed May 27 22:17:25 2015

@author: mathias
"""

import libCoal
import numpy as np
import matplotlib.pyplot as pl

class simulateKingmanFiniteSitesMutation(libCoal.simulateKingman):

    def __init__(self,n,mutationRate,L,*args):
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
        self.coal = libCoal.coalescent(libCoal.partition([[i] for i in range(n)]))
#        self.simulateCoalescent(self.args)
        self.L = L
        self.MRCA_seq = np.zeros(self.L,dtype=int)
        self.site_mutationCount = np.zeros(self.L,dtype=int)
        self.sequences_mutationCount = np.zeros(self.n, dtype=int)
        self.sequences = np.zeros((self.n,self.L),dtype=int)
        self.simulateCoalescent()

#    def preSimulationSteps(self):
#        self.L = self.args[0]
#        self.MRCA_seq = np.zeros(self.L,dtype=int)
#        self.site_mutationCount = np.zeros(self.L,dtype=int)
#        self.sequences_mutationCount = np.zeros(self.n, dtype=int)
#        self.sequences = np.zeros((self.n,self.L),dtype=int)

    def postSimulationSteps(self):
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
    def getS(self):
        return np.array(self.sequences, dtype=int)

    def getSiteMutationCounts(self):
        return np.array(self.site_mutationCount, dtype = int)

    def getSeq_mutationCount(self):
        return np.array(self.sequences_mutationCount, dtype = int )

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

    def countInconsistensies(self):
        counter = float("nan")
        affectedSites = filter(lambda i: self.site_mutationCount[i] > 0, range(self.L) )
        pass
        # TODO: Finish this
        return counter


def generate_plot_1(n,L,thetaMax,thetaMin=0,steps=20,N=100):
    h = (thetaMax - thetaMin)/steps
    thetas = np.arange(thetaMin,thetaMax,h)+h
    avgRate = np.zeros(len(thetas))
    for i,theta in enumerate(thetas):
        simulations = [simulateKingmanFiniteSitesMutation(n,float(theta)/2,L) for z in range(N)]
        rates = []
        for s in simulations:
#        for i in range(N):
#            s = simulations[i]
            minimalMutations = s.countMinimalMutations()
            actualMutations = len(s.coal.mutations)
            if minimalMutations > 0:
                rates.append( float(actualMutations) / minimalMutations )
        avgRate[i] = np.average(rates)

    label = "L = "+str(L)+" N = "+str(N)+" n ="+str(n)
    pl.plot(thetas/(2.0*L) , avgRate , color='blue' , label=label)
    pl.legend(loc='upper left')
    pl.savefig("plots/plot1__L_%i__N_%i__n_%i.pdf"%(L,N,n))

def runTests():
#    generate_plot_1(n=10,L=50,thetaMax=50,steps=50,N=1000)
    generate_plot_1(n=10,L=50,thetaMax=1.0,steps=3,N=1000)

#runTests()