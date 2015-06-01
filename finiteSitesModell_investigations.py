# -*- coding: utf-8 -*-
"""
Created on Wed May 27 22:17:25 2015

@author: mathias
"""

import libCoal
import numpy as np
import matplotlib.pyplot as pl
from scipy.special import binom

class simulator_KingmanFiniteSites(libCoal.simulateKingman):

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

    def countInconsistColumnPairs(self,verbose=False):
        counter = 0
        affectedSites = filter(lambda i: self.site_mutationCount[i] > 0, range(self.L) )
        for s1 in affectedSites:
            for s2 in filter(lambda x: x > s1 , affectedSites):
                if s1!=s2 and isInconsistentColumnPair(self.sequences[:,s1],self.sequences[:,s2]):
                    counter += 1
                    if verbose:
                        print "Sites %i and %i are inconsistent"%(s1,s2)

        return counter

def isInconsistentColumnPair(c1,c2):

    n = len(c1)

#    AX_rows = map(lambda i: c1[i] == 0 and c2[i] != 0, range(n))
#    YA_rows = map(lambda i: c1[i] != 0 and c2[i] == 0, range(n))
#    YX_rows = map(lambda i: c1[i] != 0 and c2[i] != 0, range(n))
#
#    AX_occurs = reduce(lambda x,y: x or y, AX_rows)
#    YA_occurs = reduce(lambda x,y: x or y, YA_rows)
#    YX_occurs = reduce(lambda x,y: x or y, YX_rows)

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
    pl.figure(1)
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"(actual # mutations) / (# visible mutations)")
    pl.plot(thetas/L , avgRate , color='blue' , label=label)
    pl.legend(loc='upper left')
    pl.savefig("plots/plot1__L_%i__N_%i__n_%i.pdf"%(L,N,n))

    #generate plot 2
    pl.figure(2)
#    label = "(x;y) = (theta/L ; fraction of inconsistent columns)\nL,N,n = %i,%i,%i"%(L,N,n)
    label = "L,N,n = %i,%i,%i"%(L,N,n)
    pl.xlabel(r"$\frac{\theta}{L}$")
    pl.ylabel(r"#inconsistent column-pairs / $\binom{L}{2}$")
    pl.plot(thetas/L, inconsistencyCount, color = "red", label = label)
    pl.legend(loc='upper left')
    pl.savefig("plots/plot2__L_%i__N_%i__n_%i.pdf"%(L,N,n))

    #generate plot 3
    pl.figure(3)
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

#runTests()