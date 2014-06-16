#test.py
#Copyright (C) 2014  Mathias Christensen Cronjaeger KitKat2.0@gmail.com

#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

import time
import recursionEquation as re
import libCoal as lc
from scipy.special import binom
import numpy as np

n = 15
alpha = 1.7
coalType = "xi_beta"
args = (alpha,)

print "Solving p-recursion for a Xi-beta coalescent with n=%i, alpha=%f..."%(n,round(alpha,3))
t_start = time.time()
p,g = re.p_and_g(n,'xi_beta',(alpha,))
t_stop = time.time()
print "...done! \ncalculations took %f seconds"%(round(t_stop-t_start,3))


#def jumpProb(part,n,q):
#    '''
#    calculates
#    P(initial jump is to a specific state with block-sizes given
#    by "part"); denoted p_lambda in my text.
#    "part" is here a partition of n encoded as a multiset,
#    i.e. part[i] == #i-blocks of part, and
#         sum(i * part[i]) == n
#    '''
#    m = []
#    for l in [j*[i] for i,j in enumerate(part) if j!=0]:
#        m.extend(l)
#    #do this right!
#    return re.multinomial(n,m)*re.fourWay_beta_collisionRate(n,[x for  x in m if x>1],args[0]) /q[n]
##                return 1*fourWay_beta_collisionRate(n,part,args[0]) /q[n]
#    ###JUMP PROB is calculated incorrectly!
##                return fourWay_beta_collisionRate(n,[x for  x in m if x>1],args[0]) /q[n]
#
#P_mat,q_vec = re.P_and_q(n,'xi_beta',args)
#
#P = np.zeros(n+1)
#for n1 in range(1,n):
#    for p in re.partitionsMultiset_constrained(n,n1,4):
#        P[n1] += jumpProb(p,n,q_vec)
#        


#n1 = 20
#b1 = 7
#
#pList = []
#for i,p in enumerate(re.partitionsMultiset_constrained(n,n1,4)):
#    probDist = []
#    for s in re.subpartitionsMultiset(p,b1):
#        prob = np.prod([binom(p[i],s[0][i]) for i in range(len(s[0]))])/binom(n1,b1)
#        probDist += [(s[0],prob)]
#    pList.append([(sum([x[1] for x in probDist]),len(probDist),p)]+probDist)
#    if sum([x[1] for x in probDist]) != 1.0:
#        print sum([x[1] for x in probDist]),len(probDist),p
#
#class testClass(object):
#	
#	def __init__(self,x):
#		self.foo = x
#		print "Initialized"
#  
#
#def fib(n):
#    if n == 0 or n==1:
#        return 1
#    else:
#        return fib(n-1) + fib(n-2)
#
#def listCalc(n):
#    l = []
#    listBuild(l,n,n,0)
#    return l
#
#def listBuild(l,n,N,length):
##    print n,N,l
#    if length == N-1:
#        return l.append(N)
#    else:
#        l.append(n)
#        listBuild(l,n-1,N,length+1)