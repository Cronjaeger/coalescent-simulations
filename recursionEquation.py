# -*- coding: utf-8 -*-
"""
Created on Tue May 13 13:02:17 2014

@author: mathias
"""

## Functions for computing SFS statistics explicitly.
import numpy as np
from scipy.special import binom
from itertools import product, izip, ifilter
from copy import copy

def beta_collisionRate(n,k):
    pass

def fourWay_collisionRate(n,k,lambda_collisionRate):
    r = len(k)
    s = n - sum(k)
    pass

def P(n):
    return P_kingman(n)

def P_kingman(n):
    P = np.r_[eye(1,n) , eye(n-1,n)]
    # should i just use eye(5,k=-1) (I dont think P[0,0] is ever used)
    return P

def G(P,q):
    '''
    returns an n+1xn+1 matrix (the first two rows/coulmns are inconsequential),
    such that the following equality holds:
    
    G(P,q)[n,m] == g(n,m)
    computed using the recursion:
    n>=m>1 implies g(n,m) = Sum_{k=m}^{n-1} P_{n,k}*g(k,m)

    Inputs:
    P[i,j] = P_{i+1,j+1} (transition matrix of a markov chain)
    q[i] = 1/-q_(i+2,i+2) (vacation rate of the block-counting process)
    '''
    q = np.r_[np.array([[0],[float('inf')]]),q]
#    N = P.shape[1] + 1
    N = len(q)+1
#    G = np.eye(N).dot(q)
    G = np.eye(N)
    # G[n,m] = g(n,m)
    
    #compute G under the assumption that G[m,m] == 1 for all m
    for n in range(2,N+1):
        for m in range(n-1,1,-1):
            G[n,m] = float(P[n-1,m-1:n-1].dot(G[m:n,m]))
    
    # scale row m of G by a factor of 1/-q_(m,m)
    G = G.dot(np.diag(q))
    return G
    
def g_ratio(k,G):
    '''
    returns a 1xn+1 matrix, g[:] = (0,...,0,G(k,k)/G(n,k),...,G(n-1,k)/G(n,k))
    '''
    n = G.size[1]-1
    g = np.c_[np.zeros((k)),G[k:n,k]]
    return g * G[n,k]**-1
    
def multinomial(n,m):
    '''
    n = int
    m = list of integers summing to n
    '''
    #TODO: figure out if I should re-implment binom myself
    if sum(m) != n:
        return 0
    else:
        return reduce(lambda x,y:x*y,[binom(sum(m[i:]),m[i]) for i in range(len(m)-1)])
        
def p(N):
    '''
    returns an (N+1)x(N+1)x(N+1) array p, such that
    p[n,k,b] == p^{(n)}[k,b]
    '''
    #initialize array
    p = np.zeros((N+1,N+1,N+1))

    #initial conditions are set
    for i in range(1,N+1):
        p[i,i,1] = 1.
    
    # we now iterate over n (first axis of p), and fill out the rest of p
    for n in range(1,N+1):
        #we iterate over k
        for k in range(2,n+1):
            # n1: number of blocks/lineages after first jump
            for n1 in range(k,n):
                '''
                the following for-statement deserves some explanation.
                esentially ifilter is like filter but better for fintering over LARGE collections.
                the exporession:
                    "all(earlier >= later for earlier, later in zip(seq, seq[1:])) and sum(seq) == n"
                is true if x is non-decreasing and 
                '''
#                for part in ifilter(lambda seq: all([earlier >= later for earlier, later in izip(seq, seq[1:])]) and sum(seq) == n,product(xrange(1,n-k+2),repeat=n1)):
#                    for subpart in  
                            
                        
#def partitions(n,k):
#    '''
#    returns a list of all partitions of n into exactly k parts:
#    '''
#    for i in range(1,n+1)
def partitionTest(x,n):
    '''
    Verify if a given sequence is a partition of N, sorted in descending order
    '''
    isSorted = all([earlier >= later for earlier, later in izip(x, x[1:])])
    sumEqualsN = sum(x) == n
    return isSorted and sumEqualsN

#def partitionRecursion(n,n1,i,part):
#    if i == 1:
#        return part.append(n-sum(part))
#    else:
#        for k in range(([1]+part)[-1],int((n-sum(part))//i)+1):
#            print i,part,k
#            npart = partitionRecursion(n,n1,i-1,copy(part)+[k])
#        if i==n1:
#            return npart

def partitions(n,n1):
    'Outputs the set {x : x is a partition of n into n1 parts}, n1 != 1 is presumed to hold'
    P = []
    for i in range(1,n//n1+1):
        buildPartitions((i,),P,n,n1,1,i)
    return P

def buildPartitions(part,P,n,n1,Len,Sum):
    ''''A recursive function used to generate all partitions of n into N parts (note this implementation does not handle the case n1 == 1 correctly'''
    if Len == n1-1:
        P.append(part+(n-Sum,))
    else:
        for i in range(part[-1],(n-Sum)//(n1-Len)+1):
            buildPartitions(part+(i,),P,n,n1,Len+1,Sum+i)