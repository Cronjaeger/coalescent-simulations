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

def plotBar(SFS1,SFS2,W=0.6,col=['blue','red'],labels=['series 1' , 'series 2']):
    '''
    An auxiliary function to generate bar-plots of Spectra
    '''
    n1 = len(SFS1) +1
    n2 = len(SFS2) +1
    w = W/2.0
    pl.bar(np.arange(1,n1)-w,SFS1,width=w,color=col[0],label=labels[0])
    pl.bar(np.arange(1,n2),SFS2,width=w,color=col[1],label=labels[1])
    pl.legend(loc='upper right')

def plotBar3(SFS1,SFS2,SFS3,W=0.6,col=['cyan','orange','grey'],labels=['series 1' , 'series 2' , 'series 3']):
    '''
    An auxiliary function to generate bar-plots of Spectra
    '''
    n1 = len(SFS1) +1
    n2 = len(SFS2) +1
    n3 = len(SFS3) +1
    w = W/3.0
    pl.bar(np.arange(1,n1)-1.5*w,SFS1,width=w,color=col[0],label=labels[0])
    pl.bar(np.arange(1,n2)-0.5*w,SFS2,width=w,color=col[1],label=labels[1])
    pl.bar(np.arange(1,n3)+0.5*w,SFS3,width=w,color=col[2],label=labels[2])
    pl.legend(loc='upper right')


##TEST1: Compare expected normalized SFS with empirical average of norm. SFS
##  of simulated coalescents
n = 25
theta = 4
trials = 10000
beta = 1.0
psi = 0.5
c = 1. #other values of c are probably not that relevant.
T_max = float('inf')
#coalescentType = 'kingman'
#coalescentType = 'xi_beta'
#coalescentType = 'xi_pointMass'
#coalescentType = 'xi_lambda_beta'
#coalescentType = 'lambda_beta'
coalescentType = 'xi_ew'
#coalescentType = 'lambda_ew'

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
elif coalescentType =='lambda_ew':
    for i in range(trials):
        Pi_2.append(lc.simulateLambdaEldonWakely(n,theta/2.,T_max,psi))
elif coalescentType == 'xi_ew':
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

#del Pi_1
#del Pi_2




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
#eps = 0.000001
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



##TEST 4. Plot normalized expected SFS of different coalsescents and compare them
#n = 20
#theta = 2.
##trials = 100
#eps = 10**-3
#alpha = 2-eps
#psi = 1- eps
#c = 1. # other values of c are probably not that relevant.
#T_max = float('inf')
##coalescentTypes = ['lambda_beta','xi_beta']
##coalescentTypes = ['lambda_ew', 'xi_ew']
#coalescentTypes = ['lambda_ew', 'lambda_pointmass']
#
#print "computing expected average normalized SFS from recursion equation..."
#t1 = time.time()
##Kingman
#print "kingman vs. %s vs %s, n=%i , psi=%f"%(coalescentTypes[0],coalescentTypes[1],n,psi)
#xi_K,phi_K,p_K,g_K = re.expectedSFS(n,'kingman',theta,c)
#
###beta
##xi_lam,phi_lam,p_lam,g_lam = re.expectedSFS(n,coalescentTypes[0],theta,alpha)
##xi_xi,phi_xi,p_xi,g_xi = re.expectedSFS(n,coalescentTypes[1],theta,alpha)
#
### Eldon & Wakely
#xi_lam,phi_lam,p_lam,g_lam = re.expectedSFS(n,coalescentTypes[0],theta,c,psi)
##xi_xi,phi_xi,p_xi,g_xi = re.expectedSFS(n,coalescentTypes[1],theta,c,psi)
#
## Pointmass
##xi_lam,phi_lam,p_lam,g_lam = re.expectedSFS(n,coalescentTypes[0],theta,psi)
#xi_xi,phi_xi,p_xi,g_xi = re.expectedSFS(n,coalescentTypes[1],theta,psi)
#
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
#
#print "plotting results"
#
##plot normalized spectra
#x = np.arange(1,n)
###Beta
##label1 = 'norm SFS (%s_%s )'%(coalescentTypes[0],str(round(alpha,3)))
##label2 = 'norm SFS (%s_%s )'%(coalescentTypes[1],str(round(alpha,3)))
#
## E&W or Pointmass
#label1 = 'norm SFS (%s_%s )'%(coalescentTypes[0],str(round(psi,3)))
#label2 = 'norm SFS (%s_%s )'%(coalescentTypes[1],str(round(psi,3)))
#label3 = 'norm SFS Kingman'
#lab = [label1,label2,label3]
#
###pl.plot()
##pl.plot(x , phi_lam[1:] , color='blue' , label=label1)
##pl.plot(x , phi_xi[1:] , color='red' , label=label2)
##pl.plot(x , phi_K[1:] , color='green' , label=label3)
##pl.legend(loc='upper right')
###pl.plot()
##del x
#
##plotBar(phi_lam[1:],phi_xi[1:])
#plotBar3(phi_lam[1:],phi_xi[1:],phi_K[1:],W=0.8,labels=lab)

##Test 5: Compare empirical expected spectra of Xi- versus Lambda coalescents
#n = 10
#theta = 10.
#trials = 300
#beta = 1.2
#psi = 1
#c = 1. # other values of c are probably not that relevant.
#T_max = float('inf')
#coalescentType = 'xi_beta'
#coalescentType = 'lambda_beta'
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

##this can be sped up for larger N
#print "computing expected average normalized SFS from recursion equation..."
#t1 = time.time()
#xi_1,phi_1,p_1,g_1 = re.expectedSFS(n,'lambda_beta',theta,beta)
#xi_2,phi_2,p_2,g_2 = re.expectedSFS(n,'xi_beta',theta,beta)
#t2 = time.time()
#print "done! elapsed time = %s sec \n"%(str(round(t2-t1,3)))
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
#
##check that the theoretical expected tree-ength matches the empirical one.
##compute empirical tree-length:
#lengths = np.array([x.coal.computeTreeLength() for x in Pi_2])
#empLength = np.average(lengths)
#
##compute expected theoretical tree-lenghts
#thrLength = sum(l*g[n,l] for l in range(2,n+1))
#
##print "arerage empirical tree-length:\t%f\ntheoretical tree-length:\t%f\nRelative error:\t%f"%(empLength,thrLength,abs(empLength -thrLength)/thrLength)

##Test 6
## Run the same coalescent "xi_pointmass" for different ranges of parameters
#
##coalescentType = 'xi_pointmass'
##coalescentType = 'xi_beta'
#coalescentType = 'xi_EW'
#n,theta,c = 25, 2., 1.0
#eps = 10**-9
#alphaList = [i/20.0 for i in np.arange(0,21,4)]
#alphaList[0] += eps
#alphaList[-1] += -eps
#psiList = alphaList
##psiList = [i/20. for i in range(1,21)]
##psiList = []
#K = len(psiList)
#x = np.arange(1,n)
#xi,phi,p,g,labels = [False]*K,[False]*K,[False]*K,[False]*K,[False]*K
#for i,psi in enumerate(psiList):
##    labels[i] = "Lambda = dirac_%s"%(str(round(psi,2)))
#    labels[i] = 'Xi-EW_%s'%(str(round(psi,2)))
#    print "working on %s. case %i of %i"%(labels[i],i+1,K)
##    xi[i],phi[i],p[i],g[i] = re.expectedSFS(n,coalescentType,theta,psi)
#    xi[i],phi[i],p[i],g[i] = re.expectedSFS(n,coalescentType,theta,c,psi)
#for i in range(len(psiList)):
#    pl.plot(x , phi[i][1:] , label=labels[i])
#pl.legend(loc='upper right')

#TEST 7
# compare the results obtained from recursionequation.py with those from 
# betaXiCoalescent_SFS.py
#import betaXiCoalescent_SFS as sfsExact
#coalescentType = 'xi_beta'
#n, theta, beta = 10, 2., 1.0
#
#xi1, phi1, p1, g1 = re.expectedSFS(n, coalescentType, theta, beta)
#xi2, phi2, p2, g2 = sfsExact.expectedSFS(n, coalescentType, theta, beta)
#
#label1 = 'old implementation'
#label2 = 'new implementation'
#
#x = np.arange(1,n)
#pl.plot(x , xi1[1:] , label=label1)
#pl.plot(x , xi2[1:-1] , label=label2)
#pl.legend(loc='upper right')
#
#thrLength1 = sum(l*g1[n,l] for l in range(2,n+1))
#thrLength1_1 = sfsExact.expectedTreeLength(g1)
#thrLength2 = float(sfsExact.expectedTreeLength(g2))
#
#print thrLength1,thrLength2