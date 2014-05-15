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
    P = np.r_[np.eye(1,n) , np.eye(n-1,n)]
    # should i just use eye(5,k=-1) (I dont think P[0,0] is ever used)
    return P

def q(n):
    '''
    Returns the diagonal of the q-matrix associated with the block-counting process.
    '''
    pass

def G(P,q):
    '''
    returns an n+1xn+1 matrix (the first two rows/coulmns are inconsequential),
    such that the following equality holds:
    
    G(P,q)[n,m] == g(n,m)
    computed using the recursion:
    n>=m>1 implies g(n,m) = Sum_{k=m}^{n-1} P _{n,k}*g(k,m)

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
    
    # scale row m of G by a factor of 1/-q_(m ,m)
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
        return reduce(lambda x,y:x*y,[binom(sum(m[i:]),m[i]) for i in xrange(len(m)-1)])

def collisionRate(part,n):
    #TODO IMPLEMENT THIS
    pass

def jumpProb(part,n,q):
    m = []
    for l in [j*[i] for i,j in enumerate(part) if j!=0]:
        m.extend(l)
    return multinomial(n,m)*collisionRate(part,n)/q[n]

def p(N):
    '''
    returns an (N+1)x(N+1)x(N+1) array p, such that
    p[n,k,b] == p^{(n)}[k,b]
    '''
    #compute constants:
    P = P(N)
    q = Q(N)
    G = G(P,q)
    
    #initialize array
    p = np.zeros((N+1,N+1,N+1))

    #initial conditions are set
    for i in range(1,N+1):
        p[i,i,1] = 1.
    
    myProd = np.prod
    
    # we now iterate over n (first axis of p), and fill out the rest of p
    for n in range(1,N+1):
        #we iterate over k
        for k in range(2,n+1):
            gQuotient = g_ratio(k,G)
            # n1: number of blocks/lineages after first jump
            for n1 in range(k,n):
                quotResult = gQuotient[n1]
                # we iterate over how many blocks we take from the partition we generate
                for b1 in range(1,n1-k+2):
                    b1Result = quotResult*p[n1,k,b1]
                    for p in partitionsMultiset(n,n1):
                        pResult = b1Result*jumpProb(p,n)
                        for s in subpartitionsMultiset(p,b1):
                            p[n,k,s[1]] += pResult*myProd([binom(p[i],s[0][i]) for i in xrange(s[1]) if s[0][i] != 0])
### HERE               
#                '''
#                the following for-statement deserves some explanation.
#                esentially ifilter is like filter but better for fintering over LARGE collections.
#                the exporession:
#                    "all(earlier >= later for earlier, later in zip(seq, seq[1:])) and sum(seq) == n"
#                is true if x is non-decreasing and 
#                '''
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
    '''Outputs the set {x : x is a partition of n into n1 parts}.
    The partitions are returned as tuples in lexicographical order'''
    if n1==1:
        return set([(n,)])
    else:
        P = set()
        for i in xrange(1,n//n1+1):
            buildPartitions((i,),P,n,n1,1,i)
        return P

def buildPartitions(part,P,n,n1,Len,Sum):
    ''''A recursive function used to generate all partitions of n into N parts (note this implementation does not handle the case n1 == 1 correctly'''
    if Len == n1-1:
        P.add(part+(n-Sum,))
    else:
        for i in xrange(part[-1],(n-Sum)//(n1-Len)+1):
            buildPartitions(part+(i,),P,n,n1,Len+1,Sum+i)

def partitionsMultiset(n,n1):
    '''
    Works similar to partitions(n,n1), but the partitions returned are encoded as multisets encoded as lists; e.g. the partition p=(1,1,1,2,5) of 10 would be encoded p_mul=(0,3,1,0,0,1,0,0,0,0,0), the idea being p_mul[i] == p.count(i)
    '''
    if n1 ==1:
        return set([int(j==n) for j in xrange(n+1)])
    else:
        P = set()
        for i in xrange(1,n//n1+1):
            buildPartitionsMultiset([int(j==i) for j in xrange(n+1)],i,P,n,n1,1,i)
        return P

def buildPartitionsMultiset(part,last,P,n,n1,Len,Sum):
    if Len == n1-1:
        i = n-Sum
        P.add(tuple([part[j] + int(i==j) for j in xrange(len(part))]))
    else:
        for i in xrange(last,(n-Sum)//(n1-Len)+1):
#            npart = list(part)
#            npart[i] += 1
#            buildPartitionsMultiset(tuple(npart),i,P,n,n1,Len+1,Sum+i)
            buildPartitionsMultiset(tuple([part[j] + int(i==j) for j in xrange(len(part))]),i,P,n,n1,Len+1,Sum+i)

def subpartitionsMultiset(part,b1):
    '''
    returns all subpartitions of the partitions Part (encoded as a multiset),
    that can be generated, by taking exactly b1 blocks from "part". The
    returned partitions are encoded as multisets.
    '''
    n = len(part)
    if b1==n:
        return [part]
    else:
        subP = set()
        for i in [j for j in xrange(n) if part[j] != 0]:
            buildSubPartMulti(tuple([part[j] - int(j==i) for j in xrange(n)]),tuple([int(j==i) for j in xrange(n)]),n,b1-1,i,subP)
        return subP

def buildSubPartMulti(origPart,subPart,n,toGo,Sum,subP):
    if toGo==0:
        subP.add((subPart,Sum))
    else:
        for i in [j for j in xrange(n) if origPart[j] != 0]:
            buildSubPartMulti(tuple([origPart[j] - int(j==i) for j in xrange(n)]),tuple([subPart[j] + int(j==i) for j in xrange(n)]),n,toGo-1,Sum+i,subP)