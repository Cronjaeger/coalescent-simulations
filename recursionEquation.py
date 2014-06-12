# -*- coding: utf-8 -*-
"""
Created on Tue May 13 13:02:17 2014

@author: mathias
"""

## Functions for computing SFS statistics explicitly.
import numpy as np
from scipy.special import binom
from scipy.special import beta as Beta #lower-case beta is the beta-distribution in the numpy package
#from itertools import product, izip, ifilter
from copy import copy

def lambda_beta_collisionRate(b,k,alpha):
    if k > b or k < 2:
        return 0
    else:
        return Beta(k-alpha,b-k+alpha)/Beta(2-alpha,alpha)

def fourWay_beta_collisionRate(b,k,alpha):
    ###THIS FUNCTION DOES NOT RETURN WHAT IT OUGHT TO
    ####PERHAPS IT DOES NOW?
    '''
    compute the rate of (b;k[1],...,k[len(k)];s)-collisions
    since s = b -sum(k), it does not need to be an argument
    it is assumed that all entries in the vector k are non-zero.
    '''
    k = [x for x in k if x>1] #remove all 1 and 0 entires from k
    K = sum(k) #Total number of affected blocks
    if all([i==1 for i in k]) or K > b or K < 2 :
#        print "fourWay_beta_collisionRate(%i,%i,%f) == 0"%(b,k,alpha)
        return 0
    else:
        r = len(k)
        return sum([lambda_beta_collisionRate(b,K+l,alpha)*np.prod(range(4,4-(r+l),-1))/(4.0**(K+l)) for l in range(0,4-r+1)])


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
        return int(k==2) + c*(phi**k)*((1-phi)**(b-k))

def fourWay_ew_collisionRate(b,k,c,phi):
    '''
    compute the rate of (b;k[1],...,k[len(k)];s)-collisions
    since s = b -sum(k), it does not need to be an argument
    it is assumed that all entries in the vector k are integers greater
    than 0, and that len(k) < 4
    '''
    k = [x for x in k if x>1]
    K = sum(k)
    if all([i==1 for i in k]) or K > b or K < 2 :
        return 0
    else:
        r = len(k)
        K = sum(k)
        return sum([lambda_ew_collisionRate(b,K,c,phi)*np.prod(range(4,4-(r+l),-1))/4.0**(K+l) for l in range(0,4-r+1)])
        
#        P_k = multinomial(K,k)/(4.0**K) * multinomial(len(k)-k.count(0),[k.count(i) for i in range(1,K+1)])
#        print P_k,multinomial(K,k),(4.0**K),multinomial(len(k),[k.count(i) for i in range(1,K+1)])
        # The last factor in the above, counts the number of different ways
        # to arrange the numbers k[1],...k[len(k)]
#        return P_k * lambda_ew_collisionRate(b,K,c,phi)

def P_and_q(n,coalescentType,args):
    '''
        Returns
         P: the transition matrix of the block-counting process
         q: -1*diagonal of Q-matirx of block-counting process
        
    '''
    coalescentType = str.lower(coalescentType)
    if coalescentType=='kingman':
        return P_and_q_kingman(n)
    elif coalescentType=='lambda_beta':
        return P_and_q_lambda_beta(n,args)
    elif coalescentType=='xi_beta':
        return P_and_q_xi_beta(n,args)
    elif coalescentType=='lambda_ew':
        return P_and_q_lambda_EW(n,args)
    elif coalescentType=='xi_ew':
        return P_and_q_xi_EW(n,args)
    elif coalescentType=='xi_bottleneck':
        return P_and_q_bottleneck(n,args)

def P_and_q_kingman(n):
    P = np.r_[np.eye(1,n) , np.eye(n-1,n)]
    # should i just use eye(5,k=-1) (I dont think P[0,0] is ever used)
    q = q_kingman(n)
    return P,q

def P_and_q_lambda_beta(N,args):
    alpha = args[0]
    P = np.zeros((N+1,N+1))
    q = np.zeros(N+1)
    for n in xrange(2,N+1):
        for m in xrange(1,n):
            P[n,m] = binom(n,n-m+1) * lambda_beta_collisionRate(n,n-m+1,alpha)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q
    
def P_and_q_xi_beta(N,args):
#    print "args=",args
    alpha = args[0]
    P = np.zeros((N+1,N+1))
    q = np.zeros(N+1)
    for n in xrange(2,N+1):
        for m in xrange(1,n):
            for p in partitions_constrained(n,m,4):
                #TODO: is it appropriate to multiply with multinomial(n,p)?
                # yes it is!
                P[n,m] += multinomial(n,p) * fourWay_beta_collisionRate(n,[x for x in p if x > 1],alpha)
        q[n] = sum(P[n,:])
        P[n,:] = P[n,:]/q[n]
    return P,q

def P_and_q_lambda_EW(N,args):
    c = args[0]
    phi = args[1]
    q = np.zeros(N+1)
    P = np.zeros((N+1,N+1))
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
    for n in xrange(2,N+1):
        for m in xrange(1,n):
            for p in partitions_constrained(n,m,4):
                #TODO: is it appropriate to multiply with multinomial(n,p)?
                # yes it is
                P[n,m] += multinomial(n,p) * fourWay_ew_collisionRate(n,[k for k in p if k > 1],c,phi)
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
    #TODO: it seems more appropriate to set 0^-1 to float('inf'). Does this break anything?
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
    n>=m>1 implies g(n,m) = Sum_{k=m}^{n-1} P _{n,k}*g(k,m)

    Inputs:
    P[i,j] = P_{i+1,j+1} (transition matrix of a markov chain)
    q_diag[i] = -q_(i,i) (vacation rate of the block-counting process)
    '''
#    q = np.r_[np.array([[0],[float('inf')]]),q_diag]
    q_G = copy(q_diag)
    for i,x in enumerate(q_G):
        q_G[i] = reciprocal(x)
#    N = P.shape[1] + 1
    N_G = len(q_G)
#    G = np.eye(N).dot(q)
    G_G = np.eye(N_G)
    # G[n,m] = g(n,m)
    
    #compute G under the assumption that G[m,m] == 1 for all m
    for n in range(2,N_G):
        for m in range(n-1,1,-1):
            G_G[n,m] = float(P[n-1,m-1:n-1].dot(G_G[m:n,m]))
    
    # scale row m of G by a factor of 1/-q_(m ,m)
    G_G = G_G.dot(np.diag(q_G,0))
    return G_G
    
def g_ratio(k,G):
    '''
    returns a 1xn+1 matrix, g[:] = (0,...,0,G(k,k)/G(n,k),...,G(n,k)/G(n,k))
    '''
    n = G.shape[1]-1
    g_gRatio = np.concatenate((np.zeros(k),G[k:n+1,k]),1)
    return g_gRatio * G[n,k]**(-1)
    
def multinomial(n,m):
    '''
    n = int
    m = list of integers summing to n
    '''
    #TODO: re-implment binom myself, if the standard turns out to be too slow/inacurate
    mybinom = binom
    if sum(m) != n:
        return 0
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
        'xi_beta'   (args[0] = alpha)
        'xi_ew'     (args[0] = c, args[1]=phi)
    '''
#TODO: add support for lambda-coalescents

    #compute constants:
    P_mat,q_vec = P_and_q(N,coalescentType,args)
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
    elif coalescentType=='xi_beta' or coalescentType=='xi_ew':
        if coalescentType=='xi_beta':        
            def jumpProb(part,n,q):
                '''
                calculates
                P(initial jump is to a specific state with block-sizes given
                by "part"); denoted p_lambda in my text.
                "part" is here a partition of n encoded as a multiset,
                i.e. part[i] == #i-blocks of part, and
                     sum(i * part[i]) == n
                '''
                m = []
                for l in [j*[i] for i,j in enumerate(part) if j!=0]:
                    m.extend(l)
                #do this right!
                return multinomial(n,m)*fourWay_beta_collisionRate(n,[x for  x in m if x>1],args[0]) /q[n]
#                return 1*fourWay_beta_collisionRate(n,part,args[0]) /q[n]
                ###JUMP PROB is calculated incorrectly!
#                return fourWay_beta_collisionRate(n,[x for  x in m if x>1],args[0]) /q[n]

        elif coalescentType=='xi_ew':
            def jumpProb(part,n,q):
                '''
                similar to the xi_beta-case
                '''
                m = []
                for l in [j*[i] for i,j in enumerate(part) if j!=0]:
                    m.extend(l)
                #m is an encoding of part as a sequence.
                return multinomial(n,m)*fourWay_ew_collisionRate(n,m,args[0],args[1])/q[n]
    
        #By adding np.prod to the local scope, the number of global lookups per-
        # formed in the innermost for-loop is significantly reduced.    
        myProd = np.prod
        
#        # we now iterate over n (first axis of p), and fill out the rest of p
#        for n in range(1,N+1):
#            #we iterate over k
#            for k in range(2,n+1):
#                gQuotient = g_ratio(k,G_mat)
#                # n1: number of blocks/lineages after first jump
#                for n1 in range(k,n):
#                    quotResult = gQuotient[n1]
#                    # we iterate over how many blocks we take from the partition we generate
#                    for b1 in range(1,n1-k+2):
#                        b1Result = quotResult*p_mat[n1,k,b1]
##                        if p_mat[n1,k,b1] !=0: #for testing purposes
##                            print "p_mat[n1=%i,k=%i,b1=%i]=%f" % (n1,k,b1,round(p_mat[n1,k,b1],2))
#                        for p in partitionsMultiset_constrained(n,n1,4):
#                            pResult = b1Result*jumpProb(p,n,q_vec)
##                            if jumpProb(p,n,q_vec)!=0:
##                                print "not 0; jumpProb =%f"%(jumpProb(p,n,q_vec))
##                            if pResult!=0 or b1Result!=0:
##                                print "Not 0, pResult=%f, b1Result=%f"%(pResult,b1Result)
#                            for s in subpartitionsMultiset(p,b1):
##                                product = myProd([binom(p[i],s[0][i]) for i in xrange(s[1]+1) if s[0][i] != 0])
##                                x = pResult*product
##                                if x!=0 or product!=0:
##                                    print "(x,product,pResult)=(%f,%f,%f)"%(round(x,2),round(product,2),round(pResult,2))
##                                p_mat[n,k,s[1]] += x
#                                p_mat[n,k,s[1]] += pResult*myProd([binom(p[i],s[0][i]) for i in xrange(len(s[0])) if s[0][i] != 0])/binom(n1,b1)
##                                if p_mat[n,k,s[1]] == float('inf'):
##                                    print "wtf! p_mat[%i,%i,%i] = inf"%(n,k,s[1])
##                                    print "pResult=",pResult
##                                    print "pSample s from p =",myProd([binom(p[i],s[0][i]) for i in xrange(s[1]+1) if s[0][i] != 0])/binom(n1,s[1])
##                print "p_mat[%i,%i,:]="%(n,k)+str([p_mat[n,k,:]])


##In the following, I have restructured, so that k is the inner variable. This
# should speed things up considerably. Regrettibly, the results are shit so far
        for n in range(1,N+1):
            for n1 in range(1,n):
                for p in partitionsMultiset_constrained(n,n1,4):
                    pResult = jumpProb(p,n,q_vec)
#                    JumpProb was incorrectly calculated!
                    # I suspect the problem is that q_n is incorrectly calculated
                    # this was indeed the case the problem was the lack of a copy-statement
#                    # Test
#                    if pResult > 1.0:
#                        print "P(",p," | %i) = %f \n q_vec[%i]=%f \n"%(n,pResult,b,q_vec[n])
                    for b1 in range(1,n1):
                        b1Result = pResult*(binom(n1,b1)**-1)
                        for s in subpartitionsMultiset(p,b1):
                            b = s[1]
                            sResult = b1Result*myProd([binom(p[i],s[0][i]) for i in xrange(len(s[0])) if s[0][i] != 0])
                            #what is the appropriate range for k?
#                            for k in range(2,n+1):
                            for k in [x for x in range(2,n+1) if x <= n1 and b1 <= n1 - x +1]:
#                                new =  sResult*p_mat[n1,k,b1]*(G_mat[n1,k]/G_mat[n,k])
##                                testing s[1] seems to always hold?
#                                if new == 0:
##                                    print "\n sResult =%f \n G_mat[n1,k] =%f \n p_mat[%i,%i,%i] = %f"%(sResult,G_mat[n1,k],n1,k,b,p_mat[n1,k,b])
#                                    print "(n,k,b,n1,b1,s[1])=(%i,%i,%i,%i,%i,%i)"%(n,k,b,n1,b1,s[1])
                                p_mat[n,k,b] += sResult*p_mat[n1,k,b1]*(G_mat[n1,k]/G_mat[n,k])
        return p_mat,G_mat
    
    ###CASE: Lambda-coalescents
    elif coalescentType=='lambda_beta' or coalescentType=='lambda_ew':
        for n in range(1,N+1):
            for k in range(2,n+1):
                for b in range(1,n-k+2):
                    for n1 in range(k,n):
                        res = 0.0
                        if b > n-n1:
                            res += (b-n+n1)/float(n1)*p_mat[n1,k,b-n+n1]
                        if b < n1:
                            res += (n1-b)/float(n1)*p_mat[n1,k,b]
                        p_mat[n,k,b] += (P_mat[n,n1]*G[n1,k]/G[n,k])*res
        return p_mat,G_mat

def expectedSFS(n,coalescentType,tetha,*args):
    p_mat,G_mat = p_and_g(n,coalescentType,args)
    SFS = np.zeros(n)
    normaLizedSFS = np.zeros(n)
    for i in range(1,n):
        SFS[i] = tetha/2.0 * sum([p_mat[n,k,i]*k*G_mat[n,k] for k in range(2,n-i+2)])
        normaLizedSFS[i] = SFS[i]/sum([l*G_mat(n,l)for l in range(2,n+1)])
    return SFS,normaLizedSFS
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
        return set([(n,)])
    else:
        P = set()
        for i in xrange(1,n//n1+1):
            buildPartitions((i,),i,P,n,n1,1,i)
        return P

def buildPartitions(part,last,P,n,n1,Len,Sum):
    ''''A recursive function used to generate all partitions of n into N parts (note this implementation does not handle the case n1 == 1 correctly'''
    if Len == n1-1:
        P.add(part+(n-Sum,))
    else:
        for i in xrange(last,(n-Sum)//(n1-Len)+1):
            buildPartitions(part+(i,),i,P,n,n1,Len+1,Sum+i)

def partitionsMultiset(n,n1):
    '''
    Works similar to partitions(n,n1), but the partitions returned are encoded as multisets encoded as lists; e.g. the partition p=(1,1,1,2,5) of 10 would be encoded p_mul=(0,3,1,0,0,1,0,0,0,0,0), the idea being p_mul[i] == p.count(i)
    '''
    if n1 ==1:
        return set(tuple([int(j==n) for j in xrange(n+1)]))
    else:
        P = set()
        for i in xrange(1,n//n1+1):
            buildPartitionsMultiset(tuple([int(j==i) for j in xrange(n+1)]),i,P,n,n1,1,i)
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

def partitionsMultiset_constrained(n,n1,maxBigBlocks):
    "exactly like Partitions-multiset, only the total number of non-singleton blocks is constrained"
    if n1 ==1:
        P = set()
        P.add(tuple([int(j==n) for j in xrange(n+1)]))
        return P
    else:
        P = set()
        singletonBlocks = max(n1 - maxBigBlocks,0)
        initialPart = [0,singletonBlocks]+(n-1)*[0]
        for i in xrange(1,(n-singletonBlocks)//(n1-singletonBlocks)+1):
            buildPartitionsMultiset(tuple([initialPart[j] + int(j==i) for j in xrange(n+1)]),i,P,n,n1,1+singletonBlocks,i+singletonBlocks)
        return P

def partitions_constrained(n,n1,maxBigBlocks):
    if n1 ==1:
        P = set()
        P.add((n,))
        return P
    else:
        P = set()
        singletonBlocks = max(n1 - maxBigBlocks,0)
        initialPart = singletonBlocks*[1]
        for i in xrange(1,(n-singletonBlocks)//(n1-singletonBlocks)+1):
            buildPartitions(tuple(initialPart + [i]),i,P,n,n1,1+singletonBlocks,i+singletonBlocks)
        return P

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
        subP = set()
        subP.add((part,sum([i*j for i,j in enumerate(part)])))
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