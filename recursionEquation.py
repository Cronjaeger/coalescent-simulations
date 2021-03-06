# -*- coding: utf-8 -*-
"""
Created on Tue May 13 13:02:17 2014

@author: mathias
"""
import numpy as np

import mpmath as mp
#Multiple precision arithmetics. Availiable at http://mpmath.org/
from math import factorial
from scipy.special import binom
from scipy.special import beta as Beta #lower-case beta is the beta-distribution in the numpy package
#from itertools import product, izip, ifilter
from copy import copy

def partitionTest(x,n):
    '''
    Verify if a given sequence is a partition of N, sorted in descending order
    '''
    isSorted = all([earlier >= later for earlier, later in zip(x, x[1:])])
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
    The partitions are generated in lexicographical order, but
    returned as a set (a data-structure without ordering) for
    optimization purposes'''
    if n1==1:
#        return set([(n,)])
        return [(n,)]
    else:
#        P = set()
        P = []
        for i in xrange(1,n//n1+1):
            buildPartitions((i,),i,P,n,n1,1,i)
        return P

def buildPartitions(part,last,P,n,n1,Len,Sum):
    ''''A recursive function used to generate all partitions of n into N parts (note this implementation does not handle the case n1 == 1 correctly'''
    if Len == n1-1:
#        P.add(part+(n-Sum,))
        P.append(part+(n-Sum,))
    else:
        for i in xrange(last,(n-Sum)//(n1-Len)+1):
            buildPartitions(part+(i,),i,P,n,n1,Len+1,Sum+i)

def partitionsMultiset(n,n1):
    '''
    Works similar to partitions(n,n1), but the partitions returned are encoded as multisets encoded as lists; e.g. the partition p=(1,1,1,2,5) of 10 would be encoded p_mul=(0,3,1,0,0,1,0,0,0,0,0), the idea being p_mul[i] == p.count(i)
    '''
    if n1 ==1:
#        return set(tuple([int(j==n) for j in xrange(n+1)]))
        return [tuple([int(j==n) for j in xrange(n+1)])]
    else:
#        P = set()
        P = []
        for i in xrange(1,n//n1+1):
            buildPartitionsMultiset(tuple([int(j==i) for j in xrange(n+1)]),i,P,n,n1,1,i)
        return P

def buildPartitionsMultiset(part,last,P,n,n1,Len,Sum):
    if Len == n1-1:
        i = n-Sum
#        P.add(tuple([part[j] + int(i==j) for j in xrange(len(part))]))
        P.append(tuple([part[j] + int(i==j) for j in xrange(len(part))]))
    else:
        for i in xrange(last,(n-Sum)//(n1-Len)+1):
#            npart = list(part)
#            npart[i] += 1
#            buildPartitionsMultiset(tuple(npart),i,P,n,n1,Len+1,Sum+i)
            buildPartitionsMultiset(tuple([part[j] + int(i==j) for j in xrange(len(part))]),i,P,n,n1,Len+1,Sum+i)

def partitionsMultiset_constrained(n,n1,maxBigBlocks):
    "exactly like Partitions-multiset, only the total number of non-singleton blocks is constrained"
    if n1 ==1:
#        P = set()
#        P.add(tuple([int(j==n) for j in xrange(n+1)]))
        P = []
        P.append(tuple([int(j==n) for j in xrange(n+1)]))
        return P
    else:
#        P = set()
        P = []
        singletonBlocks = max(n1 - maxBigBlocks,0)
        initialPart = [0,singletonBlocks]+(n-1)*[0]
        for i in xrange(1,(n-singletonBlocks)//(n1-singletonBlocks)+1):
            buildPartitionsMultiset(tuple([initialPart[j] + int(j==i) for j in xrange(n+1)]),i,P,n,n1,1+singletonBlocks,i+singletonBlocks)
        return P

def partitions_constrained(n,n1,maxBigBlocks):
    if n1 ==1:
#        P = set()
        P = []
#        P.add((n,))
        P.append((n,))
        return P
    else:
#        P = set()
        P = []
        singletonBlocks = max(n1 - maxBigBlocks,0)
        initialPart = singletonBlocks*[1]
        if singletonBlocks == n1-1:
            return tuple(initialPart + [n - singletonBlocks] )
        for i in xrange(1,(n-singletonBlocks)//(n1-singletonBlocks)+1):
            buildPartitions(tuple(initialPart + [i]),i,P,n,n1,1+singletonBlocks,i+singletonBlocks)
        return P

def NEWSubpartitionsMultiset(part,b1):
#    n = len(part)
#    if b1==sum(part):
#        subP = []
#        subP.append((part,sum([i*j for i,j in enumerate(part)])))
#        return subP
#    else:
#        subP = []
#        blocksGeq = list([int(part[i] > 0)*sum(part[i:]) for i in range(n)])
#        suitableChoices = (i for i,x in enumerate(blocksGeq) if x >=b1)
##        print blocksGeq
#        for i in suitableChoices:
#            subPart = tuple([int(j==i) for j in xrange(n)])
##            newBlocksGeq = [0 for j in range(i)] + [blocksGeq[i]-(blocksGeq[i-1]+1)] + blocksGeq[i+1:]
#            blocks = tuple([int(i<=j)*part[j] - int(j==i) for j in xrange(n)])
##            print subPart,'\n',newBlocksGeq,'\n'
#            NEWbuldSubpartitionsMultiset(blocks,subPart,n,b1-1,i,subP)
    n = len(part)
    n1 = sum(part)
    subP_list = []
    if n1 <= b1:
        if n1 == b1:
            subP_list.append((part,sum([i*j for i,j in enumerate(part)])))
        return subP_list
    else:
        bmax = max([i for i in range(n) if sum(part[i:]) >= b1])
        for i in (j for j in range(1,bmax+1) if part[j] > 0):
            l = tuple([int(i==j) for j in range(n)])
            a = tuple([part[j] - int(i==j) for j in range(n)])
            bmin = i
            bmax_new = max([j for j in range(bmax,n) if sum(a[j:]) >= b1 - 1])
            NEWbuldSubpartitionsMultiset(l,a,bmin,bmax_new,b1-1,n,i,subP_list)
        return subP_list
def NEWbuldSubpartitionsMultiset(l,a,bmin,bmax,toGo,n,Sum,subP_list):
    if toGo==0:
        subP_list.append((l,Sum))
    else:
        for i in (j for j in range(bmin,bmax+1) if a[j] > 0):
#            bmin_new = i
            l_new = tuple([l[j] + int(i==j) for j in range(n)])
            a_new = tuple([a[j] - int(i==j) for j in range(n)])
            bmax_new = max([j for j in range(bmax,n) if sum(a_new[j:]) >= toGo - 1])
##            print l_new,a_new,toGo,bmin_new,bmax_new
            NEWbuldSubpartitionsMultiset(l_new, a_new, i , bmax_new, toGo-1, n, Sum+i, subP_list)
#            NEWbuldSubpartitionsMultiset(l_new,a_new, i, max([j for j in range(bmax,n) if sum(a[j:]) >= toGo - 1]), toGo-1, n, Sum+i, subP_list)

#def NEWbuldSubpartitionsMultiset(blocks,subPart,n,toGo,Sum,subP):
#    if toGo==0:
#        subP.append((subPart,Sum))
#        print '\n'
#    else:
#        blocksGeq = list([int(blocks[i] > 0)*sum(blocks[i:]) for i in range(n)])
#        suitableChoices = (i for i,x in enumerate(blocksGeq) if x >=toGo)
#        for i in suitableChoices:
#            newSubPart = tuple([subPart[j] + int(j==i) for j in range(n)])
#            newBlocks = tuple([int(i<=j)*blocks[j] - int(j==i) for j in range(n)])
##            newBlocksGeq = [0 for i in range(i)]+[blocksGeq[i]-1]+blocksGeq[i+1:]
#            print newBlocks,newSubPart,i,toGo-1
#            buildSubPartMulti(newBlocks, newSubPart,n,toGo-1,Sum+i,subP)    
#    pass
def subpartitionsMultiset(part,b1):
    '''
    returns all subpartitions of the partitions Part (encoded as a multiset),
    that can be generated, by taking exactly b1 blocks from "part". The
    returned partitions are encoded as multisets.

    Each subpartiton is encoded (s,sum), where s is a multiset-encoding of the
    sub-partition, and "sum" is the sum of the block-sizes. "sum" is passed on
    as a result, so that It does not have to be calculated separately at a
    later point in time.
    '''
    n = len(part)
    if b1==n:
#        subP = set()
        subP = []
    
#        subP.add((part,sum([i*j for i,j in enumerate(part)])))
        subP.append((part,sum([i*j for i,j in enumerate(part)])))
        return subP
    else:
#        subP = set()
        subP = []
        #TODO: at the moment this mmethod adds the same subpartitions multiple times. it is not efficient. a temporary fix has been made.
        for i in (j for j in xrange(n) if part[j] != 0):
#        for i in xrange(n):
#            if part[i] != 0:
            buildSubPartMulti(tuple([part[j] - int(j==i) for j in xrange(n)]),tuple([int(j==i) for j in xrange(n)]),n,b1-1,i,subP)
        return list(set(subP)) #a hack to remove duplicate entries. The code is not very efficient

def buildSubPartMulti(origPart,subPart,n,toGo,Sum,subP):
    if toGo==0:
#        subP.add((subPart,Sum))
        subP.append((subPart,Sum))
    else:
        for i in (j for j in xrange(n) if origPart[j] != 0):
#        for i in xrange(n):
#            if origPart[i] != 0:
            buildSubPartMulti(tuple([origPart[j] - int(j==i) for j in xrange(n)]),tuple([subPart[j] + int(j==i) for j in xrange(n)]),n,toGo-1,Sum+i,subP)

def subpartitionProb(part,subpart,n1,b1,verify=True):
    '''
    IN:
        part,n1:     "part" is a (multiset-) partition (of n into n1 parts)
        subpart,b1:  "subpart" is a (multiset-) partition (of b into b1 parts)
        verify:     Should we check if subpart is a partition of part., and
                    re-calculate b1 and b2
    OUT:
        p:          Probability that one obtains "subpart" by picking b1 blocks
                    from "part"
    '''
    #verify that part subsumes subpart

    if verify:
        n1 = sum(part)
        b1 = sum(subpart)
        if not (b1<=n1 and all([x[1] <= x[0] for x in zip(part,subpart)])):
            print "%s does not subsume %s"%(str(part),str(subpart))
            return 0.0
    return np.prod([binom(x[0],x[1]) for x in zip(part,subpart)])/binom(n1,b1)

def partitionToMultiset(part):
    '''
        input: a partition encoded as a non-ascending sequence
        output: the multiset-encoding of the input-partition
        
        example:
        part = (4,2,1,1,1) partition of 9
        partitionToMultiset(part) = (0,3,1,0,1,0,0,0,0,0)
    '''
    return tuple([part.count(i) for i in range(sum(part)+1)])

def lambda_beta_collisionRate(b,k,alpha):
    if k > b or k < 2:
        return 0
    else:
        return Beta(k-alpha,b-k+alpha)/Beta(2-alpha,alpha)

def fallingFactorial(n,k):
    return np.prod(range(n,n-k,-1))

def fourWay_beta_collisionRate(b,k,alpha):
    '''
    compute the rate of (b;k[1],...,k[len(k)];s)-collisions
    since s = b -sum(k), it does not need to be an argument
    it is assumed that all entries in the vector k are non-zero.
    '''
    k = [x for x in k if x>1] #remove all 1 and 0 entires from k
    K = sum(k) #Total number of affected blocks
#    print b,k,K
    if all([i==1 for i in k]) or K > b or K < 2 :
        return 0
    else:
        r = len(k)
        s = b-K
        l_max = min(4-r,s)
        return sum([(binom(s,l) * lambda_beta_collisionRate(b,K+l,alpha))*(fallingFactorial(4,l+r)/(4.0**(K+l))) for l in range(0,l_max+1)])
#        return sum([binom(b-K,l) * lambda_beta_collisionRate(b,K+l,alpha)*np.prod(range(4,4-(r+l),-1))/(4.0**(K+l)) for l in range(0,4-r+1)])
#        rate = 0.0
#        for l in range(0,l_max+1):
#            rate += binom(s,l)*lambda_beta_collisionRate(b,K+l,alpha)*np.prod(range(4,4-(r+l),-1))/(float(4)**(K+l))
#        return rate

#        P_k = multinomial(K,k)/(4.0**K) * multinomial(len(k)-k.count(0),[k.count(i) for i in range(1,K+1)])
#        ''' The first factor counts the ways that one can partition a K-set into subsets, with sizes given by the vector k (which has), divided by the total number ow ways one can 
#         The last factor in the above, counts the number of different ways
#         to arrange the numbers k[1],...k[len(k)]'''
#        if P_k > 1: #for testing
#        print P_k,K,k
#            print (k , [k.count(i) for i in range(1,K+1)])

#        x = P_k * lambda_beta_collisionRate(b,K,alpha)
#        if x==0.:
#            print P_k,lambda_beta_collisionRate(b,K,alpha)
#        print "fourWay_beta_collisionRate(%i,%i,%f) = %f * %f = %f "%(b,k,alpha,P_k,lambda_beta_collisionRate(b,K,alpha),x)
#        return P_k * lambda_beta_collisionRate(b,K,alpha)

def lambda_ew_collisionRate(b,k,c,phi):
    if k > b or k < 2:
        return 0
    else:
        return (2./(2. + phi*phi)) * int(k==2) + c*(phi**k)*((1-phi)**(b-k))/(2. + phi*phi)

def fourWay_ew_collisionRate(b,k,c,phi):
    '''
    compute the rate of (b;k[1],...,k[len(k)];s)-collisions
    since s = b -sum(k), it does not need to be an argument
    it is assumed that all entries in the vector k are integers greater
    than 0, and that len(k) < 4
    '''
#    k = [x for x in k if x>1]
#    K = sum(k)
#    if all([i==1 for i in k]) or K > b or K < 2 :
#        return 0
#    else:
#        r = len(k)
#        K = sum(k)
#        return sum([binom(b-K,l)*lambda_ew_collisionRate(b,K,c,phi)*np.prod(range(4,4-(r+l),-1))/4.0**(K+l) for l in range(0,4-r+1)])
#        
#        P_k = multinomial(K,k)/(4.0**K) * multinomial(len(k)-k.count(0),[k.count(i) for i in range(1,K+1)])
#        print P_k,multinomial(K,k),(4.0**K),multinomial(len(k),[k.count(i) for i in range(1,K+1)])
        # The last factor in the above, counts the number of different ways
        # to arrange the numbers k[1],...k[len(k)]
#        return P_k * lambda_ew_collisionRate(b,K,c,phi)

    k = [x for x in k if x>1] #remove all 1 and 0 entires from k
    K = sum(k) #Total number of affected blocks
#    print b,k,K
    if all([i==1 for i in k]) or K > b or K < 2 :
        return 0
    else:
        r = len(k)
        s = b-K
        l_max = min(4-r,s)
        return sum([(binom(s,l) * lambda_ew_collisionRate(b,K+l,c,phi))*(fallingFactorial(4,l+r)/(4.0**(K+l))) for l in range(0,l_max+1)])

def lambda_pointMass_collisionRate(b,k,phi):
    phi = float(phi)
    if k>b or k<2:
        return 0.0
    else:
        return (phi**(k))*((1-phi)**(b-k))

def fourWay_pointMass_collisionRate(b,k,phi):
    k = [x for x in k if x>1] #remove all 1 and 0 entires from k
    K = sum(k) #Total number of affected blocks
    if all([i==1 for i in k]) or K > b or K < 2 :
        return 0.0
    else:
        r = len(k)
        return sum([binom(b-K,l) * lambda_pointMass_collisionRate(b,K+l,phi)*np.prod(range(4,4-(r+l),-1))/(4.0**(K+l)) for l in range(0,4-r+1)])
        
def P_and_q(n,coalescentType,args):
    '''
        Returns
         P: the transition matrix of the block-counting process
         q: -1*diagonal of Q-matirx of block-counting process
        
    '''
    coalescentType = str.lower(coalescentType)
    if coalescentType=='kingman' or coalescentType=='xi_kingman':
        return P_and_q_kingman(n)
    elif coalescentType=='lambda_beta' or coalescentType=='xi_lambda_beta':
        return P_and_q_lambda_beta(n,args)
    elif coalescentType=='xi_beta':
        return P_and_q_xi_beta(n,args)
    elif coalescentType=='lambda_ew':
        return P_and_q_lambda_EW(n,args)
    elif coalescentType=='xi_ew':
        return P_and_q_xi_EW(n,args)
    elif coalescentType=='xi_bottleneck':
        return P_and_q_bottleneck(n,args)
    elif coalescentType=='lambda_pointmass':
        return P_and_q_lambda_pointMass(n,args)
    elif coalescentType=='xi_pointmass':
        return P_and_q_xi_pointMass(n,args)
    else:
        print "Unknown coalescent-type"

def P_and_q_kingman(n):
    P = np.eye(n+1,k=-1)
    P[1,0] = 0.
    P[1,1] = 1.
#    P = np.r_[np.eye(1,n) , np.eye(n-1,n)]
    # should i just use eye(5,k=-1) (I dont think P[0,0] is ever used)
    q = q_kingman(n)
    return P,q

def P_and_q_lambda_beta(N,args):
    alpha = args[0]
    P = np.zeros((N+1,N+1))
    P[1,1] = 1.
    q = np.zeros(N+1)
    for n in xrange(2,N+1):
        for m in xrange(1,n):
            P[n,m] = binom(n,n-m+1) * lambda_beta_collisionRate(n,n-m+1,alpha)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q
    
def P_and_q_xi_beta(N,args):
    alpha = args[0]
    P = np.zeros((N+1,N+1))
    P[1,1] = 1.
    q = np.zeros(N+1)

#Original implementation.
    for n in xrange(2,N+1):
        for m in xrange(1,n):
            for p in partitions_constrained(n,m,4):
#            for p in partitions(n,m):
                p_mul = partitionToMultiset(p)
                k_vec = [x for x in p if x > 1]
                factor = (multinomial(n,p)/np.prod(map(factorial, p_mul)))
#                factor = multinomial(n,k_vec+[n - sum(k)]) / np.prod(map(factorial,p_mul[2:])) #from schweinsberg 2000 equation (3). (equivalent to above formula)
                P[n,m] += factor * fourWay_beta_collisionRate(n,k_vec,alpha)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q
    #    
   
#    #new attempt to reimplement
#    # begin auxiliary functions
#    def factorPart(part):
#        return multinomial(sum(part),part)/np.prod(map(factorial, partitionToMultiset(part)))
#
#    def signaturePart(part):
#        return (sum(part),[k for k in part if k > 1],sum([int(k==1) for k in part]))
#
#    def lambdaRate(b,k):
#        return lambda_beta_collisionRate(b,k,alpha)
#
#    def collisionRate(signature):
#        b,k_vec,s = signature[0], signature[1], signature[2]
#        K,r = sum(k_vec),len(k_vec)
#        l_max = min(4-r,s)
#        return sum([binom(s,l) * lambdaRate(b,K+l) * fallingFactorial(4,r+l)/(4**K+l) for l in range(l_max+1)])
#        
#    def ratePart(part):
#        return factorPart(part) * collisionRate(signaturePart(part))
#    #end auxiliary functions
#    
#    for n in range(1,N+1):
#        parts = []
#        partsBySize = []
#        for m in range(1,n):
#            partsNew = partitions_constrained(n,m,4)
#            parts += partsNew
#            partsBySize.append(partsNew)
#        q[n] = sum(map(ratePart,parts))
#        for i,partsBySize in enumerate(partsBySize):
#            P[n,i+1] = sum(map(ratePart,partsBySize))/q[n]
#    return P,q

# attempt to implement using multiple-precision arithmetics
#    P_lambdaTest,q_lambdaTest = P_and_q_lambda_beta(N,args)
#    decimalPlaces = 100
#    with mp.workdps(decimalPlaces):
#        alpha = args[0]
#        P = np.zeros((N+1,N+1))
#        P[1,1] = 1.
#        q = np.zeros(N+1)
#        Q_mat = [ [ [] for j in range(N+1) ] for i in xrange(N+1)]
#        for n in range(2,N+1):
#            for m in range(1,n):
#                for p in partitions_constrained(n,m,4):
##                    q_p = multinomial(n,p,multiplePrecision=True) * mp.mpf(fourWay_beta_collisionRate(n,[x for x in p if x > 1],alpha))
#                    p_mul = tuple([p.count(i) for i in range(n+1)])
#                    q_p = multinomial(n,p,multiplePrecision=True) * mp.mpf(fourWay_beta_collisionRate(n,[x for x in p if x > 1],alpha)) / np.prod(map(mp.fac,p_mul))
#                    Q_mat[n][m].append(q_p)
#            q_n = sum(map(sum,Q_mat[n]))
#            q[n] = float(q_n)
#
#            if q[n] > q_lambdaTest[n]:
#                print "Error! For n=%i q_xi > q_lambda\n\tq_xi=%f\n\tq_lambda=%f\n"%(n,q[n],q_lambdaTest[n])
#            for m in range(1,n):
#                p_nm = sum(Q_mat[n][m])/(q_n)
##                p_nm = sum(Q_mat[n][m])/(q_n*mp.fac(n-2))
#                P[n,m] = float(p_nm)
    return P,q

def P_and_q_lambda_EW(N,args):
    c = args[0]
    phi = args[1]
    q = np.zeros(N+1)
    P = np.zeros((N+1,N+1))
    P[1,1] = 1.
    for n in xrange(1,N+1):
        for m in xrange(1,n):
            ### P_and_q_lambda_EW, Is q[n] correctly calculated?
            P[n,m] = binom(n,n-m+1)*lambda_ew_collisionRate(n,n-m+1,c,phi)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q    

def P_and_q_xi_EW(N,args):
    c = args[0]
    phi = args[1]
    q = np.zeros(N+1)
    P = np.zeros((N+1,N+1))
    P[1,1] = 1.
#    for n in xrange(2,N+1):
#        for m in xrange(1,n):
#            for p in partitions_constrained(n,m,4):
#                P[n,m] += multinomial(n,p) * fourWay_ew_collisionRate(n,[k for k in p if k > 1],c,phi)
#        q[n] = sum(P[n,:])
#        P[n,:] = P[n,:]/q[n]
#    return P,q
    for n in xrange(2,N+1):
        for m in xrange(1,n):
            for p in partitions_constrained(n,m,4):
                p_mul = partitionToMultiset(p)
                k_vec = [x for x in p if x > 1]
                factor = (multinomial(n,p)/np.prod(map(factorial, p_mul)))
                P[n,m] += factor * fourWay_ew_collisionRate(n,k_vec,c,phi)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q

def P_and_q_lambda_pointMass(N,args):
    phi = args[0]
    P = np.zeros((N+1,N+1))
    P[1,1] = 1.
    q = np.zeros(N+1)
    for n in xrange(2,N+1):
        for m in xrange(1,n):
            P[n,m] = binom(n,n-m+1) * lambda_pointMass_collisionRate(n,n-m+1,phi)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q

def P_and_q_xi_pointMass(N,args):
    phi = args[0]
    q = np.zeros(N+1)
    P = np.zeros((N+1,N+1))
    P[1,1] = 1.
#    for n in xrange(2,N+1):
#        for m in xrange(1,n):
#            for p in partitions_constrained(n,m,4):
#                P[n,m] += multinomial(n,p) * (np.prod(map(factorial,partitionToMultiset(p)))**-1) * fourWay_pointMass_collisionRate(n,[k for k in p if k > 1],phi)
#        q[n] = sum(P[n,:])
#        P[n,:] = P[n,:]/q[n]

    for n in xrange(2,N+1):
        for m in xrange(1,n):
            for p in partitions_constrained(n,m,4):
                p_mul = partitionToMultiset(p)
                k_vec = [x for x in p if x > 1]
                factor = (multinomial(n,p)/np.prod(map(factorial, p_mul)))
                P[n,m] += factor * fourWay_pointMass_collisionRate(n,k_vec,phi)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q

def P_and_q_bottleneck(N,args):
    #TODO: implement this
    pass

def q(n,coalescentType,args=[]):
    '''
    Returns -q_(i,i), where q is the Q-matrix associated with the block-counting process of the coalescent started from n blocks.
    q[0] and q[1] are both set to 0 (they should play no role in smulations)
    '''
    coalescentType = str.lower(coalescentType)
    if coalescentType=='kingman':
        return q_kingman(n)
#    elif coalescentType=='lambda_beta':
#        return q_lambda_beta(n,args)
#    elif coalescentType=='xi_beta':
#        return q_xi_beta(n,args)
#    elif coalescentType=='lambda_ew':
#        return q_lambda_EW(n,args)
#    elif coalescentType=='xi_ew':
#        return q_xi_EW(n,args)
#    elif coalescentType=='xi_bottleneck':
#        return q_bottleneck(n,args)
def q_kingman(n):
    q = np.zeros(n+1)
    for i in xrange(2,n+1):
        q[i] = binom(i,2)
    return q

#def P_and_q_lambda_beta(n,args):
#    alpha = args[0] 
#    q = np.zeros(n+1)
#    for b in xrange(2,n+1):
#        Sum = 0
#        for k in xrange(2,b):
#            Sum += binom(b,k)*lambda_beta_collisionRate(b,k,alpha)
#        q[b] = Sum
#    return q

#def q_xi_beta(n,args):
#    pass

def reciprocal(x):
    '''
    IN : x (number)
    OUT: x_inv, where x_inv = x^-1 if x!=0; x_inv = 0 if x=0
    '''
#    TODO: it seems more appropriate to set 0^-1 to float('inf'). Does this break anything?
    if x==0:
#        return float('inf')
        return 0
    else:
        return x**-1

def G(P,q_diag):
    '''
    returns an n+1xn+1 matrix (the first two rows/coulmns are inconsequential),
    such that the following equality holds:
    
    G(P,q_diag)[n,m] == g(n,m)
    computed using the recursion:
    g(m,m) = 1/-q_{m,m}
    n>=m>1 implies g(n,m) = Sum_{k=m}^{n-1} P_{n,k}*g(k,m)

    Inputs:
    P[i,j] = P_{i,j} (transition matrix of a markov chain)
    q_diag[i] = -q_(i,i) (vacation rate of the block-counting process)
    '''

    N_G = len(q_diag)

    # Comupte diagonal elements of g, using G[n,n] = 1/abs(q_n,n,)
    q_G = copy(q_diag)
    for i,x in enumerate(q_G):
        q_G[i] = reciprocal(x)
    G_G = np.diag(q_G)

    for n in range(2,N_G):
#        for m in range(n-1,1,-1):
        for m in range(2,n):
            G_G[n,m] = float(P[n,m:n].dot(G_G[m:n,m]))
#            G_G[n,m] = sum([P[n,l]*G_G[l,m] for l in range(m,n)])
    
#     scale row m of G by a factor of 1/-q_(m ,m)
#    G_G = G_G.dot(np.diag(q_G,0))
    return G_G

#    G = np.zeros((N_G,N_G))
#    for n,x in enumerate(q_diag):
#        if x!= 0:
#            G[n,n] = float(x)**-1
#
#    for n in range(2,N_G):
#        for m in range(2,n):
#            for k in range(m,n):
#                G[n,m] += P[n,k]*G[k,m]
#    
#    return G
    
def g_ratio(k,G):
    '''
    returns a 1xn+1 matrix, g[:] = (0,...,0,G(k,k)/G(n,k),...,G(n,k)/G(n,k))
    '''
    n = G.shape[1]-1
    g_gRatio = np.concatenate((np.zeros(k),G[k:n+1,k]),1)
    return g_gRatio * G[n,k]**(-1)
    
def multinomial(n,m,multiplePrecision=False,decimalPlaces=40):
    '''
    n = int
    m = list of integers summing to n
    '''
    mybinom = binom
    if multiplePrecision: #use multiple-precision arithmetics
        with mp.workdps(decimalPlaces):
            if sum(m) != n:
                return mp.mpf('0.0')
            else:
                return np.prod([mp.mpf(mybinom(sum(m[i:]),m[i])) for i in xrange(len(m)-1)])
    else: #use floating-point arithmetics
        if sum(m) != n:
            return 0.0
        else:
    #        if len(m)==1: m = list(m)+[0] #else the reduce statement does not work
    #        return reduce(lambda x,y:x*y,[mybinom(sum(m[i:]),m[i]) for i in xrange(len(m)-1)])
            return np.prod([mybinom(sum(m[i:]),m[i]) for i in xrange(len(m)-1)])
    


##OLD
#
#def collisionRate(part,n):
#    pass
#

def p_and_g(N,coalescentType,args):
    '''
    returns an (N+1)x(N+1)x(N+1) array p, such that
    p[n,k,b] == p^{(n)}[k,b]
    supported coalescent types are:
        'kingman'
        'xi_beta'           (args[0] = alpha)
        'xi_ew'             (args[0] = c, args[1]=phi)
        'xi_pointMass'      (args[0] = phi)
        'lambda_beta'       (args[0] = alpha)
        'lambda_ew'         (args[0] = c, args[1]=phi)
        'lambda_pointMass'  (args[0] = phi)
        
    '''

    #compute constants:
    P_mat,q_vec = P_and_q(N,coalescentType,args)
#    print "1: q_%s=%s"%(coalescentType , str([round(q,3) for q in q_vec]))
    G_mat = G(P_mat,q_vec)
    #initialize array
    p_mat = np.zeros((N+1,N+1,N+1))

    #initial conditions are set
    for i in range(1,N+1):
        p_mat[i,i,1] = 1.

    coalescentType = str.lower(coalescentType)

    if coalescentType=='kingman':
        for n in range(1,N+1):
            for k in range(2,n+1):
                for b in range(1,n-k+2):
                    p_mat[n,k,b] = binom(n-b-1,k-2)/binom(n-1,k-1)
        return p_mat,G_mat
    
    #the case of four-way Xi-coalescents is treated
    elif coalescentType in set(('xi_beta','xi_ew','xi_pointmass','xi_kingman','xi_lambda_beta')):
        if coalescentType=='xi_beta':
            jumpProb = jumpProb_xiBet
        elif coalescentType=='xi_ew':
            jumpProb = jumpProb_xiEW
        elif coalescentType=='xi_pointmass':
            jumpProb = jumpProb_xiPointMass
        elif coalescentType=='xi_kingman':
            jumpProb = jumpProb_xiKingman
        elif coalescentType=='xi_lambda_beta':
            jumpProb = jumpProb_xi_lambda_beta
            args = (P_mat,)
    
        #By adding functions to the local scope, the number of global lookups per-
        # formed in the innermost for-loop is significantly reduced.    
        myProd = np.prod
#        mySubpartitionProb = subpartitionProb
        
#        # we now iterate over n (first axis of p), and fill out the rest of p
#        for n in range(1,N+1):
#            #we iterate over k
#            for k in range(2,n+1):
##                gQuotient = g_ratio(k,G_mat)
#                # n1: number of blocks/lineages after first jump
#                for n1 in range(k,n):
#                    quotResult = G_mat[n1,k]/G_mat[n,k]
##                    quotResult = gQuotient[n1]
#                    # we iterate over how many blocks we take from the partition we generate
#                    for b1 in range(1,n1-k+2):
#                        b1Result = quotResult*p_mat[n1,k,b1]
#                        for p in partitionsMultiset_constrained(n,n1,4):
#                            pResult = b1Result*jumpProb(p,n,q_vec)
#                            for s in subpartitionsMultiset(p,b1):
#                                p_mat[n,k,s[1]] += pResult*myProd([binom(p[i],s[0][i]) for i in xrange(len(s[0])) if s[0][i] != 0])/binom(n1,b1)

#In the following, I have restructured, so that k is the inner variable. This
# should speed things up considerably. Regrettibly, this has not halped improve performance this far.
        for n in range(1,N+1):
            for n1 in range(1,n):
                for p in partitions_constrained(n,n1,4):
                    p_mul = partitionToMultiset(p)
                    pResult = jumpProb(p,p_mul,n,q_vec,args)
                    for b1 in range(1,n1):
                        b1Result = pResult/binom(n1,b1)
#                        kRange = [x for x in range(2,n+1) if x <= n1 and b1 <= n1 - x + 1]
#                        for subPart,b in subpartitionsMultiset(p_mul,b1):
                        for subPart,b in NEWSubpartitionsMultiset(p_mul,b1):
#                            b = s[1]
                            sResult = b1Result*myProd([binom(p_mul[i],subPart[i]) for i in xrange(len(subPart)) if subPart[i] != 0])
#                            sResult = pResult*mp.mpf(mySubpartitionProb(p_mul,subPart[0],n1,b1,verify=True))
#                            sResult = pResult*mySubpartitionProb(p_mul,subPart,n1,b1,verify=True)
#                            for k in kRange:
                            for k in range(2,n+1):
#                                p_mat[n,k,s[1]] += sResult*mp.mpf(p_mat[n1,k,b1])*(mp.mpf(G_mat[n1,k])/mp.mpf(G_mat[n,k]))
                                p_mat[n,k,b] += sResult*p_mat[n1,k,b1]*G_mat[n1,k]/G_mat[n,k]

## AN ATEMPT AT IMPLEMENTING THE RECURSION FORMULA NAIVELY (in order to check
## for errors) Should be at least an order of magnitude slower than above
## implementations.
#        for n in range(1,N+1):
#            for k in range(2,n+1):
#                for b in range(1,n-k+2):
#                    for n1 in range(k,n):
#                        n1Res = G_mat[n1,k]/G_mat[n,k]
#
#                        for b1 in range(1,min(b,n1-k+1)+1):
#                            b1Res = n1Res*p_mat[n1,k,b1]/binom(n1,b1)
#
#                            for lam in partitions_constrained(n,n1,4):
#                                lam_multi = partitionToMultiset(lam)
#                                lamRes = b1Res*jumpProb(lam,lam_multi,n,q_vec,args)
###                                ##testing
###                                jProb = jumpProb(lam,n,q_vec,args)
###                                if jProb == 0.0:
###                                    print "jumpProb = %f\n\tlam_multiset=%s\n\tn,q[n]=%i,%f"%(jProb,str(lam),n,q_vec[n])
#                                for lam1 in [x for x in subpartitionsMultiset(lam_multi,b1) if x[1]==b]:
#                                    p_mat[n,k,b] += lamRes*myProd([binom(lam_multi[i],lam1[0][i]) for i in xrange(len(lam1[0])) if lam1[0][i] != 0])
        return p_mat,G_mat
    
    ###CASE: Lambda-coalescents
    elif coalescentType in set(('lambda_beta','lambda_ew','lambda_pointmass')):
        for n in range(1,N+1):
            for k in range(2,n+1):
                for b in range(1,n-k+2):
                    for n1 in range(k,n):
                        res = 0.0
                        if b > n-n1:
                            res += (b-n+n1)/float(n1)*p_mat[n1,k,b-n+n1]
                        if b < n1:
                            res += (n1-b)/float(n1)*p_mat[n1,k,b]
                        p_mat[n,k,b] += (P_mat[n,n1]*G_mat[n1,k]/G_mat[n,k])*res
        return p_mat,G_mat

def jumpProb_xiBet(part,partMul,n,q_vec,args):
    '''
    calculates
    P(initial jump is to a state with block-sizes given
    by "part"); denoted p_lambda in my text.
    "part" is here a partition of n encoded as a multiset,
    i.e. part[i] == #i-blocks of part, and
         sum(i * part[i]) == n
    '''
#    m = []
#    for l in [j*[i] for i,j in enumerate(part) if j!=0]:
#        m.extend(l)

#Works in floating-point arithmetics
#    return (multinomial(n,partMul)/q_vec[n])*fourWay_beta_collisionRate(n,[x for  x in part if x>1],args[0])
    return (multinomial(n,part)*fourWay_beta_collisionRate(n,[x for  x in part if x>1],args[0]))/(np.prod(map(factorial,partMul))*q_vec[n])

#    modified to work with multiple Precision arithmetics
#    return (multinomial(n,part,multiplePrecision=True)*mp.mpf(fourWay_beta_collisionRate(n,[x for  x in part if x>1],args[0])))/(np.prod(map(mp.fac,partMul))*q_vec[n])

def jumpProb_xiEW(part,partMul,n,q_vec,args):
    '''
    similar to the xi_beta-case
    '''
#    m = []
#    for l in [j*[i] for i,j in enumerate(part) if j!=0]:
#        m.extend(l)
#    #m is an encoding of part as a sequence.
#    return multinomial(n,m)*fourWay_ew_collisionRate(n,m,args[0],args[1])/q_vec[n]
    return (multinomial(n,part)*fourWay_ew_collisionRate(n,[x for  x in part if x>1],args[0],args[1]))/(np.prod(map(factorial,partMul))*q_vec[n])
    
def jumpProb_xiPointMass(part,partMul,n,q_vec,args):
    '''
    similar to the xi_beta-case
    '''
#    m = []
#    for l in [j*[i] for i,j in enumerate(part) if j!=0]:
#        m.extend(l)
#    #m is an encoding of part as a sequence.
#    return multinomial(n,p)*fourWay_pointMass_collisionRate(n,p,args[0])/(q_vec[n] * np.prod(map(factorial,p_mul)))
    return (multinomial(n,part)*fourWay_pointMass_collisionRate(n,[x for  x in part if x>1],args[0]))/(np.prod(map(factorial,partMul))*q_vec[n])

def jumpProb_xiKingman(part,n,q_vec,args):
    '''
        calculates jump-probabilities of kingmans coalescent seen as a
        xi -coalescent
    '''
    if part[1]==n-2 and part[2]==1:
        return 1.
    else:
        return 0

def jumpProb_xi_lambda_beta(p,p_mul,n,q_vec,args):
    '''
        Calculates the probability that the first jump of a beta coalescent is
        to the partition "part".
        used to check recursion-formula.
        args[0] = P
    '''
    if sum(p_mul[2:])==1:
        k = p_mul[2:].index(1) + 2
        n1 = n - k + 1
        return args[0][n,n1]
    else:
        return 0

def jumpProbTest(n,coalescentType,args,outputDist=False):
    '''
    verify that sum_(lam \in {partitions of n}) P(first jump from n to lam) ==1
    Used to guage the effect of rounding errors and testing.
    '''

    P,q_vec = P_and_q(n,coalescentType,args)

    if coalescentType=='xi_beta':
        jumpProb = jumpProb_xiBet
    elif coalescentType=='xi_ew':
        jumpProb = jumpProb_xiEW
    elif coalescentType=='xi_pointmass':
        jumpProb = jumpProb_xiPointMass
    elif coalescentType=="xi_kingman":
        jumpProb = jumpProb_xiKingman
    elif coalescentType == 'xi_lambda_beta':
        jumpProb = jumpProb_xi_lambda_beta
        args = (P,)
        
    l = []
    for n1 in range(1,n):
        for lam in list(partitions_constrained(n,n1,4)):
            l.append((lam,float(jumpProb(lam,partitionToMultiset(lam),n,q_vec,args))))
#            print "woo",n,lam,jumpProb(lam,n,q_vec,args)
    if outputDist:
        return l
    else:
        s = [x[1] for x in l]
        return sum(s),np.average(s),[x for x in l if x[1]==min(s) or x[1] ==max(s)]

def expectedSFS(n,coalescentType,tetha,*args):
    '''The function to be called from outside this program. It returns the
    following four arrays:
    - The expected site-frequency-spectrum
    - the expected normalized site-frequency spectrum
    - the solution of the p(n)[k,b]-recursion equations
    - the solution of the g(n,k)-recursions
    '''
    p_mat,G_mat = p_and_g(n,coalescentType,args)
    SFS = np.zeros(n)
    normaLizedSFS = np.zeros(n)
    normFactor = sum([l*G_mat[n,l] for l in range(2,n+1)])*tetha/2.0
    for i in range(1,n):
        SFS[i] = tetha/2.0 * sum([p_mat[n,k,i]*k*G_mat[n,k] for k in range(2,n-i+2)])        
        normaLizedSFS[i] = SFS[i]/normFactor
    return SFS,normaLizedSFS,p_mat,G_mat