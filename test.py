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

import recursionEquation as re
from scipy.special import binom
import numpy as np

n = 50
n1 = 20
b1 = 7

pList = []
for i,p in enumerate(re.partitionsMultiset_constrained(n,n1,4)):
    probDist = []
    for s in re.subpartitionsMultiset(p,b1):
        prob = np.prod([binom(p[i],s[0][i]) for i in range(len(s[0]))])/binom(n1,b1)
        probDist += [(s[0],prob)]
    pList.append([(sum([x[1] for x in probDist]),len(probDist),p)]+probDist)
    if sum([x[1] for x in probDist]) != 1.0:
        print sum([x[1] for x in probDist]),len(probDist),p

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