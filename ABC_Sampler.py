# -*- coding: utf-8 -*-
"""
ABC_Sampleer.py
 Implement a basic ABC-campler for computing posterior distributions on
 parameters of interest of a coalescent modell.

Created on Mon Jan 26 15:46:11 2015

@author: mathias
"""

import libCoal as lc #Personal library for simulating coalescents
import recursionEquation as re # personal script from MSc-Thesis
import numpy as np
import matplotlib.pyplot as pl
import time

"""
Begin: Auxiliary functions and variables
"""

def L2(obsData,simData,computeForObs=False):

    if computeForObs:
        obsSFS = obsData.coal.computeNormalizedSFS()
    else:
        obsSFS = obsData

    simSFS = simData.coal.computeNormalizedSFS()
    return np.sqrt(sum(np.square(obsSFS - simSFS)))

def L1(obsData,simData,computeForObs=False):
    if computeForObs:
        obsSFS = obsData.coal.computeNormalizedSFS()
    else:
        obsSFS = obsData

    simSFS = simData.coal.computeNormalizedSFS()
    return sum(map(abs,obsSFS - simSFS))

def L2_sorted(obsData,simData,computeForObs=False):

    if computeForObs:
        obsSFS = obsData.coal.computeNormalizedSFS()
    else:
        obsSFS = obsData

    simSFS = simData.coal.computeNormalizedSFS()
    return np.sqrt(sum(np.square(np.sort(obsSFS) - np.sort(simSFS))))

def diffTreeLength(obsData,simData,computeForObs=False):
    return abs(obsData.coal.computeTreeLength() - simData.coal.computeTreeLength())

def sampleBeta(a,b):
    return float(np.random.beta(a,b,1))

def priorUniform1to2():
    """Returns a floating point number uniformly distributed on [1,2)"""
    return float(1+np.random.random_sample())

def priorDiscrete(choices = [1.0, 2 - 2**-20],weights = [0.5,0.5]):
    return np.random.choice(a=choices,p=weights)
    

class ABC_Sampler(object):
    """
    Implements a basic ABC-sampler.
    """
    
    def __init__(self,data,mutationRate = 1.0):
        """
        A New Sampler is initialized using 
        """
        self.data = data
        self.mutationRate = mutationRate
#        self.metric = metericpass
#        self.statistic = statistic
#        self.prior = prior
    
    def samplePosteriorBeta_new(self,a,b):
        pass
    
    def samplePosteriorBeta(self,N=1,pseudometric=L2,prior=priorUniform1to2,epsilon=0.1,n_leaves = 100):
        
        obsData = self.data
        alphas = []
        acceptedData = []
        simulationCounter = [0]

        def test(simData):
            return pseudometric(obsData,simData) <= epsilon


        def sample(counter = simulationCounter):
            """Sample new Parameter"""
            alphaNew = prior() # sample from prior

            """Sample synthetic Data (Beta-coalescent)"""
            P,q = re.P_and_q_lambda_beta(n_leaves,(alphaNew,))
            coalescentNew = lc.simulateLambdaBeta(n_leaves,self.mutationRate,float('inf'),alphaNew,P,q)

            counter[0] += 1 #increment simulation counter            
            
            return alphaNew,coalescentNew
        
        """Take N approximate postrior-samples"""
        for i in xrange(N):
            alphaNew, coalescentNew = sample()
            
            while not test(coalescentNew):
                alphaNew, coalescentNew = sample()
            
            alphas.append(alphaNew)
            acceptedData.append(coalescentNew)
        
        acceptanceRate = N/float(simulationCounter[0])
        
        return alphas,acceptedData,acceptanceRate

def testSampler(N = 10**2, n_leaves = 100, mutationRate = 5.0, epsilon = 1.0, prior=priorUniform1to2, pseudometric=L2, obsData = False):

    if not obsData:
        underlyingCoalescent = lc.simulateKingman(n_leaves,mutationRate*50,float('inf')) #(boost mutationRate)
        obsData = underlyingCoalescent.coal.computeNormalizedSFS()

    mySampler = ABC_Sampler(obsData,mutationRate)
    
    print "Taking %s posterior samples with epsilon = %s"%(str(N),str(round(epsilon,3)))
    t1 = time.time()
    alphas,data,acceptanceRate = mySampler.samplePosteriorBeta(N,pseudometric,prior,epsilon,n_leaves)
    t2 = time.time()
    
    print "elapsed time = %s sec \n"%(str(round(t2-t1,3)))
    pl.hist(alphas,bins=min(50,max(N/10,1)))
    
    return alphas,data,acceptanceRate