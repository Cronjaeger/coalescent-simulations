#! /usr/bin/python2.7
# -*- coding: utf-8 -*-

"""
Created on Wed May  7 18:11:33 2014

@author: mathias
"""
import time
import recursionEquation as re
import libCoal as lc
import numpy as np
import matplotlib.pyplot as pl

#Parameters for simulation
#n = 50
#theta = 1
#trials = 1000
#phi = 0.08
#T_max = float('inf')
#
#Pi_EW = []
#Pi_EW4 = []
#
#for i in range(trials):
#    Pi_EW.append(lc.simulateLambdaEldonWakely(n,theta,T_max,phi))
#    Pi_EW4.append(lc.simulateLambdaEldonWakely_FourWay(n,theta,T_max,phi))
#
#EW_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_EW]) for i in range(n)])
#EW4_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_EW4]) for i in range(n)])
#
#x = np.arange(1,n+1)
#
#pl.plot(x,EW_SFS_AVG, color='blue',label='SFS of EW where theta=1, phi=0.8')
#pl.plot(x , EW4_SFS_AVG,color='red',label='SFS of 4-way EW where theta=1, phi=0.8')
#pl.legend(loc='upper right')
#

#n = 50
#theta = 2.
#trials = 1000
#beta = 1.5
#T_max = float('inf')
#
#Pi_1 = []
#Pi_2 = []
#
#for i in range(trials):
#    Pi_1.append(lc.simulateLambdaBeta(n,theta,T_max,beta))
#    Pi_2.append(lc.simulateLambdaBeta_FourWay(n,theta,T_max,beta))
#
##Bet_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_1]) for i in range(n)])
##Bet4_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_2]) for i in range(n)])
#
#Pi1_normSFS = [x.coal.computeNormalizedSFS() for x in Pi_1]
#Pi2_normSFS = [x.coal.computeNormalizedSFS() for x in Pi_2]
#
#Pi1_normSFS_AVG = np.array([np.average([x[i] for x in Pi1_normSFS]) for i in range(n)])
#Pi2_normSFS_AVG = np.array([np.average([x[i] for x in Pi2_normSFS]) for i in range(n)])
#
#x = np.arange(1,n+1)
#label1 = 'b-coal. ; theta=1, beta=1.5'
#label2 = '4-way b-coal. ; theta=1, beta=1.5'
#pl.plot(x , map(np.log,Pi1_normSFS_AVG) , color='blue' , label=label1)
#pl.plot(x , map(np.log,Pi2_normSFS_AVG) , color='red' , label=label2)
#pl.legend(loc='upper right')
#
##pl.plot(x,[np.log(y) for y in Bet_SFS_AVG], color='blue',label='b-coal. ; theta=1, beta=1.5 (log-scale)')
##pl.plot(x ,[np.log(y) for y in Bet4_SFS_AVG],color='red',label='4-way b-coal. ; theta=1, beta=1.5 (log-scale)')
#pl.legend(loc='upper right')
#
#del Pi_1
#del Pi_2

n = 20
theta = 1.
trials = 2000
beta = 1.5
T_max = float('inf')

Pi_1 = []
Pi_2 = []

print "Running %i simulations of a xi_beta coalescent, with n=%i, theta=%s, alpha=%s..."%(trials,n,str(round(theta,3)),str(round(beta,3)))
t1 = time.time()
for i in range(trials):
#    Pi_1.append(lc.simulateLambdaBeta(n,theta,T_max,beta))
    Pi_2.append(lc.simulateLambdaBeta_FourWay(n,theta,T_max,beta))
t2 = time.time()
print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#Bet_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_1]) for i in range(n)])
Bet4_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_2]) for i in range(n)])

#Pi1_normSFS = [x.coal.computeNormalizedSFS() for x in Pi_1]

print "computing average normalized SFS..."
t1= time.time()
Pi2_normSFS = [x.coal.computeNormalizedSFS() for x in Pi_2]

#Pi1_normSFS_AVG = np.array([np.average([x[i] for x in Pi1_normSFS]) for i in range(n)])
Pi2_normSFS_AVG = np.array([np.average([x[i] for x in Pi2_normSFS]) for i in range(n)])
t2 = time.time()
print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))

print "computing expected average normalized SFS from recursion equation..."
t1 = time.time()
xi,phi,p,g = re.expectedSFS(n,"xi_beta",theta,beta)
t2 = time.time()
print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))

#print "plotting results"
x = np.arange(1,n+1)
##plot spectra
#label1 = 'expected SFS (recursions)'
#label2 = 'emperical mean of SFS'
#pl.plot(x , np.array(list(xi[1:])+[0.]) , color='blue' , label=label1)
#pl.plot(x , Bet4_SFS_AVG , color='red' , label=label2)

##plot normalized spectra
#label1 = 'expected norm. SFS (recursions)'
#label2 = 'emperical mean of norm. SFS'
#pl.plot(x , np.array(list(phi[1:])+[0.]) , color='blue' , label=label1)
#pl.plot(x , Pi2_normSFS_AVG , color='red' , label=label2)
#pl.legend(loc='upper right')
del x

#check that the theoretical expected tree-ength matches the empirical one.
#compute empirical tree-length:
lengths = np.array([x.coal.computeTreeLength() for x in Pi_2])
empLength = np.average(lengths)

#compute expected theoretical tree-lenghts
thrLength = sum(l*g[n,l] for l in range(2,n+1))

print "arerage empirical tree-length:\t%f\ntheoretical tree-length:\t%f\nRelative error:\t%f"%(empLength,thrLength,(empLength -thrLength)/thrLength)

#del Pi_1
#del Pi_2