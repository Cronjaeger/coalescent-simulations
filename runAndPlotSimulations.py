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

#TEST1: Compare expected normalized SFS with empirical average of norm. SFS
#  of simulated coalescents
n = 5
theta = 4.
trials = 300
beta = 1.0
psi = .3
c = 1. # other values of c are probably not that relevant.
T_max = float('inf')
coalescentType = 'xi_beta'
#coalescentType = 'xi_pointMass'
#coalescentType = 'xi_lambda_beta'

Pi_1 = []
Pi_2 = []

coalescentType = coalescentType.lower()

print "Running %i simulations of a %s coalescent, with n=%i, theta=%s, alpha=%s, phi=%s..."%(trials,coalescentType,n,str(round(theta,3)),str(round(beta,3)),str(round(psi,3)))
t1 = time.time()
if coalescentType == 'xi_beta':
    P,q = re.P_and_q_lambda_beta(n,(beta,))
    for i in range(trials):
        Pi_2.append(lc.simulateLambdaBeta_FourWay(n,theta/2.,T_max,beta,P,q))
elif coalescentType == 'lambda_beta' or coalescentType=='xi_lambda_beta':
    P,q = re.P_and_q_lambda_beta(n,(beta,))
    for i in range(trials):
        Pi_2.append(lc.simulateLambdaBeta(n,theta/2.,T_max,beta,P,q))
elif coalescentType =='xi_ew':
    for i in range(trials):
        Pi_2.append(lc.simulateLambdaEldonWakely(n,theta/2.,T_max,psi))
elif coalescentType == 'lambda_ew':
    for i in range(trials):
        Pi_2.append(lc.simulateLambdaEldonWakely_FourWay(n,theta/2.,T_max,psi))
elif coalescentType == 'kingman' or coalescentType == 'xi_kingman':
    for i in range(trials):
        Pi_2.append(lc.simulateKingman(n,theta/2.,T_max))
elif coalescentType == 'lambda_pointmass':
    for i in range(trials):
        Pi_2.append(lc.simulateLambdaPoint(n,theta/2.,T_max,psi))
elif coalescentType == 'xi_pointmass':
    for i in range(trials):
        Pi_2.append(lc.simulateLambdaPoint_FourWay(n,theta/2.,T_max,psi))
t2 = time.time()
print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))

print "computing expected average normalized SFS from recursion equation..."
t1 = time.time()
if coalescentType in set(['xi_beta','lambda_beta','xi_lambda_beta']):
    xi,phi,p,g = re.expectedSFS(n,coalescentType,theta,beta)
elif coalescentType == 'xi_ew' or coalescentType == 'lambda_ew':
    xi,phi,p,g = re.expectedSFS(n,coalescentType,theta,c,psi)
elif coalescentType == 'kingman' or coalescentType == 'xi_kingman':
    xi,phi,p,g = re.expectedSFS(n,coalescentType,theta)
elif coalescentType == 'xi_pointmass' or coalescentType == 'lambda_pointmass':
    xi,phi,p,g = re.expectedSFS(n,coalescentType,theta,psi)
t2 = time.time()
print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))

print "computing average normalized SFS..."
t1= time.time()
#Bet_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_1]) for i in range(n)])
Bet4_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_2]) for i in range(n)])

factor = sum([l*g[n,l] for l in range(2,n+1)])*theta/2.0
#Pi1_normSFS = [x.coal.computeNormalizedSFS() for x in Pi_1]
Pi2_normSFS = [x.coal.computeNormalizedSFS(factor) for x in Pi_2]

#Pi1_normSFS_AVG = np.array([np.average([x[i] for x in Pi1_normSFS]) for i in range(n)])
Pi2_normSFS_AVG = np.array([np.average([x[i] for x in Pi2_normSFS]) for i in range(n)])
t2 = time.time()
print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))


print "plotting results"
x = np.arange(1,n)

#plot normalized spectra
label1 = 'expected norm. SFS (recursions)'
label2 = 'emperical mean of norm. SFS'
pl.plot(x , phi[1:] , color='blue' , label=label1)
pl.plot(x , Pi2_normSFS_AVG[:-1] , color='red' , label=label2)
#pl.plot(x , xi[1:] , color='blue' , label=label1)
#pl.plot(x , Bet4_SFS_AVG[:-1] , color='red' , label=label2)
pl.legend(loc='upper right')
del x

#check that the theoretical expected tree-ength matches the empirical one.
#compute empirical tree-length:
lengths = np.array([x.coal.computeTreeLength() for x in Pi_2])
empLength = np.average(lengths)

#compute expected theoretical tree-lenghts
thrLength = sum(l*g[n,l] for l in range(2,n+1))

print "arerage empirical tree-length:\t%f\ntheoretical tree-length:\t%f\nRelative error:\t%f"%(empLength,thrLength,abs(empLength -thrLength)/thrLength)
#
##del Pi_1
##del Pi_2

##TEST 2:
## simulate and run kingman-coalescents
#n = 15
#theta = 2
#trials = 2000
#T_max = float('inf')
#coalescentType = 'kingman'
#
#print "Running %i simulations of a %s coalescent, with n=%i, theta=%s ..."%(trials,coalescentType,n,str(round(theta,3)))
#Pi_1 = []
#Pi_2 = []
#t1 = time.time()
#for i in range(trials):
#    Pi_2.append(lc.simulateKingman(n,theta/2,T_max,1.))
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#print "computing expected average normalized SFS from recursion equation..."
#t1 = time.time()
#xi,phi,p,g = re.expectedSFS(n,coalescentType,theta)
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#print "computing emperical average normalized SFS..."
#t1= time.time()
#Bet4_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_2]) for i in range(n)])
#
#factor = sum([l*g[n,l] for l in range(2,n+1)])*theta/2.0
#Pi2_normSFS = [x.coal.computeNormalizedSFS(factor) for x in Pi_2]
#
#Pi2_normSFS_AVG = np.array([np.average([x[i] for x in Pi2_normSFS]) for i in range(n)])
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#print "plotting results"
#x = np.arange(1,n)
#
##plot normalized spectra
#label1 = 'expected norm. SFS (recursions)'
#label2 = 'emperical mean of norm. SFS'
#pl.plot(x , phi[1:] , color='blue' , label=label1)
#pl.plot(x , np.array(Pi2_normSFS_AVG[:-1]) , color='red' , label=label2)
#pl.legend(loc='upper right')
#del x
#
##check that the theoretical expected tree-ength matches the empirical one.
##compute empirical tree-length:
#lengths = np.array([x.coal.computeTreeLength() for x in Pi_2])
#empLength = np.average(lengths)
#
##compute expected theoretical tree-lenghts
#thrLength = sum(l*g[n,l] for l in range(2,n+1))
#
#print "arerage empirical tree-length:\t%f\ntheoretical tree-length:\t%f\nRelative error:\t%f"%(empLength,thrLength,abs(empLength -thrLength)/thrLength)

##TEST 3: verify that lambda-beta is similar to kingman when alpha ~2
#n = 10
#eps = 0.001
#theta = 2.0
#alpha = 2.0 - eps
#
#print "computing expected average normalized SFS from recursion equation..."
#t1 = time.time()
##xi_bet,phi_bet,p_bet,g_bet = re.expectedSFS(n,'lambda_beta',theta,alpha)
#xi_bet,phi_bet,p_bet,g_bet = re.expectedSFS(n,'xi_beta',theta,alpha)
#xi_K,phi_K,p_K,g_K = re.expectedSFS(n,'kingman',theta)
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#print "plotting results"
#x = np.arange(1,n)
##plot normalized spectra
#label1 = 'norm. SFS (xi_beta_%s )'%(str(round(alpha,3)))
#label2 = 'norm. SFS (Kingman)'
#pl.plot(x , phi_bet[1:] , color='blue' , label=label1)
#pl.plot(x , phi_K[1:] , color='red' , label=label2)
#pl.legend(loc='upper right')
#del x

##TEST 4. Plot normalized expected SFS of a xi-versus-lambda coalescent
#n = 10
#theta = 1.
#trials = 100
#alpha = 1.8
#psi = 1
#c = 1. # other values of c are probably not that relevant.
#T_max = float('inf')
#
#print "computing expected average normalized SFS from recursion equation..."
#t1 = time.time()
#xi_lam,phi_lam,p_lam,g_lam = re.expectedSFS(n,'lambda_beta',theta,alpha)
#xi_xi,phi_xi,p_xi,g_xi = re.expectedSFS(n,'xi_beta',theta,alpha)
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#print "plotting results"
#x = np.arange(1,n)
##plot normalized spectra
#label1 = 'norm SFS (lambda_beta_%s )'%(str(round(alpha,3)))
#label2 = 'norm SFS (xi_beta_%s )'%(str(round(alpha,3)))
#pl.plot(x , phi_lam[1:] , color='blue' , label=label1)
#pl.plot(x , phi_xi[1:] , color='red' , label=label2)
#pl.legend(loc='upper right')
#del x


##Test 5: Compare expected spectra of Xi- versus Lambda coalescents
#n = 30
#theta = 10.
#trials = 300
#beta = 1.2
#psi = 1
#c = 1. # other values of c are probably not that relevant.
#T_max = float('inf')
##coalescentType = 'xi_beta'
##coalescentType = 'lambda_beta'
#
#Pi_1 = []
#Pi_2 = []
#
#print "Running %i simulations of a %s coalescent, with n=%i, theta=%s, alpha=%s..."%(trials,'Beta-coalescents',n,str(round(theta,3)),str(round(beta,3)))
#t1 = time.time()
#P,q = re.P_and_q_lambda_beta(n,(beta,))
#for i in range(trials):
#    Pi_1.append(lc.simulateLambdaBeta(n,theta/2.,T_max,beta,P,q))
#    Pi_2.append(lc.simulateLambdaBeta_FourWay(n,theta/2.,T_max,beta,P,q))
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
###this can be sped up for larger N
##print "computing expected average normalized SFS from recursion equation..."
##t1 = time.time()
##xi_1,phi_1,p_1,g_1 = re.expectedSFS(n,'lambda_beta',theta,beta)
##xi_2,phi_2,p_2,g_2 = re.expectedSFS(n,'xi_beta',theta,beta)
##t2 = time.time()
##print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#print "computing average normalized SFS..."
#t1= time.time()
#Bet_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_1]) for i in range(n)])
#Bet4_SFS_AVG = np.array([np.average([x.SFS[i] for x in Pi_2]) for i in range(n)])
#
##factor_1 = sum([l*g_1[np_xiBe,l] for l in range(2,n+1)])*theta/2.0
##factor_2 = sum([l*g_2[n,l] for l in range(2,n+1)])*theta/2.0
#factor_1, factor_2 = False, False
#Pi1_normSFS = [x.coal.computeNormalizedSFS(factor_1) for x in Pi_1]
#Pi2_normSFS = [x.coal.computeNormalizedSFS(factor_2) for x in Pi_2]
#
#Pi1_normSFS_AVG = np.array([np.average([x[i] for x in Pi1_normSFS]) for i in range(n)])
#Pi2_normSFS_AVG = np.array([np.average([x[i] for x in Pi2_normSFS]) for i in range(n)])
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#
#print "plotting results"
#x = np.arange(1,n)
#
##plot normalized spectra
#label1 = 'norm. SFS (Lambda_%s)'%(str(round(beta,3)))
#label2 = 'norm. SFS (4-way-Lambda_%s)'%(str(round(beta,3)))
#pl.plot(x , Pi1_normSFS_AVG[:-1] , 'b-', label=label1)
#pl.plot(x , Pi2_normSFS_AVG[:-1] , 'r-', label=label2)
#pl.legend(loc='upper right')
#del x

#check that the theoretical expected tree-ength matches the empirical one.
##compute empirical tree-length:
#lengths = np.array([x.coal.computeTreeLength() for x in Pi_2])
#empLength = np.average(lengths)
#
##compute expected theoretical tree-lenghts
#thrLength = sum(l*g[n,l] for l in range(2,n+1))

#print "arerage empirical tree-length:\t%f\ntheoretical tree-length:\t%f\nRelative error:\t%f"%(empLength,thrLength,abs(empLength -thrLength)/thrLength)

#Test 6
# Run the same coalescent "xi_pointmass" for different ranges of parameters
#
#coalescentType = 'xi_pointmass'
#n,theta = 10,2.
#psiList = [i/20. for i in range(1,21)]
##psiList = []
#K = len(psiList)
#x = np.arange(1,n)
#xi,phi,p,g,labels = [False]*K,[False]*K,[False]*K,[False]*K,[False]*K
#for i,psi in enumerate(psiList):
#    xi[i],phi[i],p[i],g[i] = re.expectedSFS(n,coalescentType,theta,psi)
#    labels[i] = "Lambda = dirac_%s"%(str(round(psi,2)))
#for i in range(len(psiList)):
#    pl.plot(x , phi[i][1:] , label=labels[i])


    
    