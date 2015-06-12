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

    def untillFirstTwoInconsistencies(self,computeS_seq = False):

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
        while inconsistencyCount < 2 and mutationCounter < M:

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

def simulateUntillTwoMutations(N = 1000, n = 20, L = 100, mutRate = 200,printFirst10 = False):
    K_list = [simulator_KingmanFiniteSites(n,mutRate,L,False) for i in xrange(N)]
    totalTypeCount = np.array((0,0,0,0))
    misses = 0
    k_res_List = []

    for K in K_list:
        res = K.untillFirstTwoInconsistencies()
        if res["Inconsistencies"] == 2:
            totalTypeCount += res["typeCount_arr"]
            k_res_List.append([K,res])
        else:
            misses += 1

    if printFirst10:
        for K in K_list[:10]:
            print chronology(K)
            print ""

    return totalTypeCount,N-misses,k_res_List

def chronology(K):
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
    typeCounts,N_eff,k_res_list = simulateUntillTwoMutations(N = N, n = n , L = L, mutRate=theta, printFirst10= False)
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


def generatePlot_of_mutationTypes(N = 1000,L = 100, n = 20, printFirsrst10 = False,show = False):

    theta = 1.2 * L
    typeCounts,N_eff,K_list = simulateUntillTwoMutations(N = N, n = n , L = L, mutRate=theta, printFirst10= printFirsrst10)

    #run simulations
    typeCounts,N_eff = simulateUntillTwoMutations(N = N, n = n , L = L, mutRate=theta, printFirst10= printFirsrst10)[:2]

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
    pl.draw()
    filename_str = "plots/bars_stoppedProcess/bar_typeFrequencies_N_%i_L_%i_n_%i"%(N_eff,L,n)
    try:
        pl.savefig(filename_str+".pdf")
        pl.savefig(filename_str+".png")
        pl.savefig(filename_str+".ps")
        pl.savefig(filename_str+".svg")
        pl.savefig(filename_str+".eps")
    except Exception:
        print "could not save in all formats"
    if show:
        pl.show()


def run_generatePlot_of_mutationTypes(arglist = [(1000,100,20)]):
    for args in arglist:
        N,L,n = args
        print "Generating plots (barcharts) for N,L,n = %i,%i,%i"%tuple(args)
        generatePlot_of_mutationTypes(N,L,n,False)

def run_generateScatterplots(arglist = [(100,100,20)]):
    for args in arglist:
        N,L,n = args
        print "Generating plots (scatterplts) for N,L,n = %i,%i,%i"%tuple(args)
        scatterplot_index_and_time_of_abnormal_mutations(N,L,n)