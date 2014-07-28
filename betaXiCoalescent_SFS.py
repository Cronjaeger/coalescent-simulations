# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 11:29:01 2014

@author: mathias

This is a seccond attempt at implementing the SFS-formula for xi-coalescents
arising from Beta coalescents split four ways.

Everything is implemented using multiple precision arithmetics
"""

import numpy as np
import mpmath as mp
import recursionEquation as re
#from scipy.special import binom
#from scipy.special import beta as Beta

alphaDefault = 1.0
mp.dps = 20 #sets number of significant figures (base 10) to be used in computations

def betaMergerRate(b, k, alpha=alphaDefault):
    if 1 < k and k <= b:
        return mp.beta(k-1, b-k+1) / mp.beta(2-alpha, alpha)
    else:
        print "wrong parameters in betaMergerRate. b, k = %s, %s"%(str(b), str(k))
        return mp.mpf('0')

def fourWayMergerRate(b, k_vec, s, lambdaMergerRate=betaMergerRate):
    K = sum([k for k in k_vec if k > 1])
    r = len([k for k in k_vec if k > 1])
    if K + s != b or K<2 or r>4:
#        print "invalid signature in fourWayMergerRate. (b, k_vec, s) = %s"%(str((b, k_vec, s)))
        return mp.mpf(0)
    else:
        l_max = min(s, 4-r)
        lRange = range(0, l_max+1)
        def summand(l):
            return mp.binomial(s, l) * ( fallingFactorial(4, r+l) / mp.power(4, K+l) ) * lambdaMergerRate(b, K+l)
        return sum([summand(l) for l in lRange])

def Q_matrix_Xi(transitionRatesList):
    '''
    compute the Q-matrix of the block-counting-process of a XI-coalescent
    the input isintended to be a list of the form:
    L[i] = transitionRatesNumberPartitions(i, ....)
    '''
    n = len(transitionRatesList)-1
    Q = np.zeros((n+1, n+1), dtype=mp.mpf)
    

    '''
    We first compute the off-diagonal elements. of Q.
    since i<j => Q_ij = 0
    '''
    for i, ratesFrom_i in enumerate(transitionRatesList[1:], start=1): #iterate over rows of Q
        for j, ratesFrom_i_to_j in enumerate(ratesFrom_i[1:], start=1): #iterate over j < i
            Q[i, j] = sum([x[2] for x in ratesFrom_i_to_j])
    
    '''
    We now compute the diagonal elements of Q:
    '''
    for i in range(1, n+1):
        Q[i, i] = -sum(Q[i, :i])
    
    return Q

def P_matrix(Q):
    P = np.zeros(Q.shape, dtype=mp.mpf)
    P[1, 1] = mp.mpf('1')
    for i in range(2, P.shape[0]):
        P[i, :i] = Q[i, :i]/(-Q[i, i])
    return P

def g_matrix(P, Q):
    g = np.zeros(Q.shape, dtype=mp.mpf)
    N = Q.shape[0] - 1 
    for n in range(2, N+1):
        g[n, n] = -1/Q[n, n]
        for m in range(2, n):
            g[n, m] = sum([P[n, k]*g[k, m] for k in range(m, n)])
    return g

def p_recursions(N, coalescentType, args):
    if coalescentType == 'xi_beta':
        alpha = args[0]
        def lambdaMergerRate(b, k):
            return betaMergerRate(b, k, alpha)
        def fourwayMergerRate(b, k_vec, s):
            return  fourWayMergerRate(b, k_vec, s, lambdaMergerRate)

    L = [transitionRatesNumberPartitions(n, fourwayMergerRate, 4) for n in range(1, N+1)]
    Q = Q_matrix_Xi([[]]+L)
    P = P_matrix(Q)
    g = g_matrix(P, Q)

    p = np.zeros((N+1, N+1, N+1), dtype=mp.mpf)
    p[1, 1, 1] = mp.mpf('1')
    
#    print L, "\n"
    for n, jumpsFrom_n in enumerate(L[1:], start=2):
        p[n, n, 1] = mp.mpf('1')
        for n1, jumpsFrom_n_to_n1 in enumerate(jumpsFrom_n[1:-1], start=1):
#            print n, n1, jumpsFrom_n_to_n1
            for lam_multi, rate in [( re.partitionToMultiset(x[0]),  x[2] ) for x in jumpsFrom_n_to_n1]:
                jumpProb = -rate/Q[n, n]
                for b1 in range(1, n1):
                    for lamSub, b in re.subpartitionsMultiset(lam_multi, b1):
                        lamSubFactor = re.subpartitionProb(lam_multi, lamSub, n1, b1)
                        for k in range(2, min(n1-b1+1, n-b+1) +1):
                            kFactor = p[n1, k, b1] * g[n1, k]/g[n, k]
                            p[n, k, b] += jumpProb * kFactor * lamSubFactor
                            print n,n1,k,b,b1,'\n',lam_multi,'\n',lamSub,'\n'
#        print '\n'
    #for testing if things add up
#    for n, l1 in enumerate(p):
#        for k, l2 in enumerate(l1):
#            if sum(l2)>0:
#                print (n, k, sum(l2))

    return p, g, Q, P,L
    
def expectedTreeLength(g):
    n = g.shape[0]-1
    return sum([l*g[n, l] for l in range(2, n+1)])

def SFS(p, g,  theta=mp.mpf('2')):
    n = g.shape[0]-1
    xi = np.array([mp.mpf('0')]+[theta/mp.mpf('2') * sum([p[n, k, i]*k*g[n, k] for k in range(2, n-i+2)]) for i in range(1, n+1)])
    return xi

def normSFS(p, g):
    theta=mp.mpf('2')
    xi = SFS(p, g, theta)
    treeLength = expectedTreeLength(g)
    psi = xi/(theta/mp.mpf('2') * treeLength)
    return psi

def expectedSFS(N, coalescentType, theta, *args):
    '''
    this is a re-implementation of the method expectedSFS from
    recursionEquation.py
    '''
    p, g, Q, P = p_recursions(N, coalescentType, args)
    xi = SFS(p,g,theta)
    phi = normSFS(p,g)
    return map(toFloat,(xi, phi, p, g))


'''
# Begin auxiliary functions #
'''
def transitionRatesNumberPartitions(n, mergerRatesCoalescent, maxNoBigblocks=4):
    '''
        returns a list lists of tuples (part, n1, rate), where:
        output[i] contains all tuples of the form (.., i, ..) 
        (Note that this implies output[0] = [])
            part - a partition of n. All partitions of n are in the list.
            n1 - the number of blocks in part.
            rate - the rate of jumps of the form:
                    from <1^n> = (1, ..., 1) ; partition of n into sigleton-blocks
                    to "part" ; some other partition of n
    '''
    
    def rate(partition):

        b, k_vec, s = computeSignature(partition)
        partition_multiset = re.partitionToMultiset(partition)

        mergerRate = mergerRatesCoalescent(b, k_vec, s)
        factorSizes = np.prod(map(mp.fac, partition))
        factorCounts = np.prod(map(mp.fac, partition_multiset))

#        print factorSizes, factorCounts, mergerRate

        return (mp.fac(n) / factorSizes*factorCounts ) * mergerRate
        
    partitionsByBlockcount = [re.partitions_constrained(n, n1, maxNoBigblocks) for n1 in range(1, n+1)]

    output = [[] for i in range(n+1)]
    for n1, partitionList in enumerate(partitionsByBlockcount, start=1):
        l = [(p, n1, rate(p)) for p in partitionList]        
        output[n1].extend(l)
            
    return output

def transitionRatesList(n, mergerRatesCoalescent=fourWayMergerRate):
    return [[]]+[transitionRatesNumberPartitions(i, mergerRatesCoalescent) for i in range(1, n+1)]

def fallingFactorial(n, k):
    return np.prod(map(mp.mpf, range(n, n-k, -1)))

def computeSignature(partition):
    b = sum(partition)
    s = partition.count(1)
    k_vec = [k for k in partition if k > 1]
    if b != s+sum(k_vec):
        print "computeSignature returns an invalid signature"
    return b, k_vec, s

def toFloat(x):
    return np.array(x,dtype='float')

'''
# End auxiliary functions #
'''