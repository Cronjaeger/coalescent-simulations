# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 08:36:25 2015

@author: mathias
"""

from finiteSitesModell_investigations import simulator_KingmanFiniteSites
# import finiteSitesModell_investigations as fsm
import sys
import numpy as np
from psiFromCSV import psiFromCSV

def simData(N,theta,L):

    simTree = simulator_KingmanFiniteSites(N,theta/2,L)

    S_redundantRowsAndColumns = simTree.getS()

    #remove null-rows
    S1 = withoutNullColumns(S_redundantRowsAndColumns)

    #remove redundant rows, and associate each row with a count-vector n instead
    S,n = S_and_n(S1)

    return S,n

def simData_XExtraMutations(N,L,k,mutationsMustBeVisible =True):

    simTree = simulator_KingmanFiniteSites(N,10*L/2.0,L,False)
    result = simTree.untillFirstXInconsistencies(X = k)

    if not mutationsMustBeVisible:
        # assure that we have k more mutations on the underlying tree than
        # we have segregating sites
        while result["Inconsistencies"] < k:
            result = simTree.untillFirstXInconsistencies(X = k)
    else:
        # Assure that we have k visible inconsistencies with the infinite sites
        # model
        while result["typeCount"][0] + result["typeCount"][2] < k:
            result = simTree.untillFirstXInconsistencies(X = k)

    S_redundantRowsAndColumns = simTree.getS()

    #remove null-rows
    S1 = withoutNullColumns(S_redundantRowsAndColumns)

    #remove redundant rows, and associate each row with a count-vector n instead
    S,n = S_and_n(S1)

    return S,n

def simData_k_mutations_total(N,L,k):
    simTree = simulator_KingmanFiniteSites(N,10.0*k,L,False)
    result = simTree.until_k_mutations(k)

    #re-run simulations if we didn't get enough mutations.
    while k > len(result['coalescent'].mutations):
        simTree = simulator_KingmanFiniteSites(N,10.0*k,L,False)
        result = simTree.until_k_mutations(k)

    S_redundant_rows_and_columns = result['S']
    S_redundant_columns, Nr = S_and_n( S_redundant_rows_and_columns )
    S_transpose, Nc = S_and_n( np.transpose(S_redundant_columns) )
    S = np.transpose(S_transpose)
    return S,Nr,Nc


def withoutNullColumns(S):
    rows,columns = S.shape
    nonNullColumns = [i for i in range(columns) if max(S[:,i]) > 0]
    columns_new = len(nonNullColumns)

    if columns_new > 0:
        S_new = S[:,nonNullColumns]
        S_new.shape = (rows,columns_new)
    else:
        S_new = np.zeros((rows,1))

    return S_new


def S_and_n(S):

    rows,columns = S.shape
    row_counts = {}

    for row in S:
        rowAsTuple = tuple(row)
        if rowAsTuple in row_counts:
            row_counts[rowAsTuple] += 1
        else:
            row_counts[rowAsTuple] = 1

    rowList = row_counts.keys()
    rowCounts = row_counts.values()

    S_new = np.array(np.r_["0,2",rowList],dtype = int)
    n_vec = np.array(rowCounts,dtype = int)

    return S_new,n_vec

def toCSV(S,Nr,Nc,fileName = False):
    '''
    Prints the configuration (S,Nr,Nc) to fileName as a .csv-file.
    If no filename is provided, prints to std-out.
    '''
    a = np.r_[np.c_[np.matrix(Nc), 0] , np.c_[S,Nr]]

    #print a

    # #sort rows so that Nr is in non-ascending order
    # a = a[a[:,-1].argsort()[::-1]]
    #
    # #sort columns such that Nc is in non-ascending order
    # a = a[:,a[0,:].argsort()[::-1]]

    if fileName == False:
        print '\n'.join([', '.join([str(a[i,j]) for j in xrange(a.shape[1])]) for i in xrange(a.shape[0])])
    else:
        np.savetxt(fileName, a, fmt = '%d', delimiter=", ")


def toCSV_old(S,n,fileName):
    a = np.c_[S,n]

    #sort rows so that n is in non-ascending order
    a = a[a[:,-1].argsort()[::-1]]

    np.savetxt(fileName, a, fmt = '%d', delimiter=", ")

# def psiFromCSV(fileName):
#     """
#     Takes the path of a .csv-file as input. The last column is presumed to
#     encode the "n"-vector enumerating the occurrance of different haplotypes;
#     The other columns are presuemd to encode haplotypes.
#
#     If the below were the context of myPhi.csv:
#     "
#     0, 1, 1, 0, 0, 2, 0, 0, 3, 3, 22
#     0, 0, 2, 1, 0, 3, 0, 0, 2, 1, 18
#     1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 12
#     0, 1, 1, 0, 2, 2, 0, 0, 3, 3, 6
#     0, 1, 1, 0, 0, 2, 0, 0, 0, 3, 6
#     1, 2, 1, 0, 0, 0, 2, 0, 0, 0, 4
#     1, 1, 1, 3, 0, 3, 0, 0, 3, 3, 3
#     1, 1, 1, 0, 0, 2, 0, 0, 3, 3, 3
#     0, 1, 1, 2, 0, 2, 0, 0, 3, 3, 3
#     1, 0, 1, 1, 2, 2, 0, 0, 0, 0, 3
#     0, 1, 1, 0, 2, 2, 0, 1, 3, 3, 3
#     0, 3, 1, 0, 0, 2, 0, 0, 3, 3, 3
#     3, 1, 1, 0, 0, 2, 0, 2, 3, 3, 3
#     0, 3, 1, 1, 0, 3, 0, 0, 0, 1, 2
#     1, 0, 1, 0, 0, 2, 1, 0, 0, 0, 2
#     0, 0, 1, 1, 2, 3, 0, 0, 2, 1, 2
#     0, 0, 1, 1, 0, 3, 0, 0, 0, 1, 1
#     0, 0, 2, 1, 0, 3, 0, 0, 0, 1, 1
#     1, 0, 1, 1, 2, 0, 0, 0, 0, 0, 1
#     1, 1, 1, 3, 0, 3, 0, 3, 3, 3, 1
#     3, 2, 2, 1, 0, 3, 0, 0, 0, 1, 1
#     "
#     Then psiFromCSV(fileName) would return S,n whereby:
#
#     >>> S
#     array([[0, 1, 1, 0, 0, 2, 0, 0, 3, 3],
#        [0, 0, 2, 1, 0, 3, 0, 0, 2, 1],
#        [1, 0, 1, 0, 0, 2, 0, 0, 0, 0],
#        [0, 1, 1, 0, 2, 2, 0, 0, 3, 3],
#        [0, 1, 1, 0, 0, 2, 0, 0, 0, 3],
#        [1, 2, 1, 0, 0, 0, 2, 0, 0, 0],
#        [1, 1, 1, 3, 0, 3, 0, 0, 3, 3],
#        [1, 1, 1, 0, 0, 2, 0, 0, 3, 3],
#        [0, 1, 1, 2, 0, 2, 0, 0, 3, 3],
#        [1, 0, 1, 1, 2, 2, 0, 0, 0, 0],
#        [0, 1, 1, 0, 2, 2, 0, 1, 3, 3],
#        [0, 3, 1, 0, 0, 2, 0, 0, 3, 3],
#        [3, 1, 1, 0, 0, 2, 0, 2, 3, 3],
#        [0, 3, 1, 1, 0, 3, 0, 0, 0, 1],
#        [1, 0, 1, 0, 0, 2, 1, 0, 0, 0],
#        [0, 0, 1, 1, 2, 3, 0, 0, 2, 1],
#        [0, 0, 1, 1, 0, 3, 0, 0, 0, 1],
#        [0, 0, 2, 1, 0, 3, 0, 0, 0, 1],
#        [1, 0, 1, 1, 2, 0, 0, 0, 0, 0],
#        [1, 1, 1, 3, 0, 3, 0, 3, 3, 3],
#        [3, 2, 2, 1, 0, 3, 0, 0, 0, 1]])
#     >>> n
#     array([22, 18, 12,  6,  6,  4,  3,  3,  3,  3,  3,  3,  3,  2,  2,  2,  1,
#         1,  1,  1,  1])
#     """
#     raw = np.genfromtxt(fileName, delimiter=',')
#     n = np.array(raw[:,-1], dtype = int)
#     S = np.array(raw[:,:-1], dtype = int)
#     return S,n

def main_old():
    #Read from input
    N = int(sys.argv[1])
    theta = float(sys.argv[2])
    L = int(sys.argv[3])
    if len(sys.argv) > 4:
        save = True
        path = sys.argv[4]
        thetaStr = ("%1.4f"%theta).replace('.','pt')
        fileName = path+"/psi__N_%i_theta_%s_L_%i.csv"%(N,thetaStr,L)
    else:
        save = False

    # simTree = finiteSitesModell_investigations.simulator_KingmanFiniteSites(N,theta/2,L)
    #
    # S_redundantRowsAndColumns = simTree.getS()
    #
    # S1 = withoutNullColumns(S_redundantRowsAndColumns)
    #
    # S,n = S_and_n(S1)

    S,n = simData(N,theta,L)

    if save:
        toCSV(S,n,fileName)
        print "Output saved to %s"%fileName

        # ##Test
        # Snew,nnew = psiFromCSV(fileName)
        # print "successfully loaded from %s"%fileName
        # print "S =\n",Snew,"\nn=\n",nnew

    else:
        print "S =\n",S,"\nn=\n",n

def main():
    #Read from input
    N = int(sys.argv[1])
    L = int(sys.argv[2])
    k = int(sys.argv[3])
    if len(sys.argv) > 4:
        save = True
        path = sys.argv[4]
        #thetaStr = ("%1.4f"%theta).replace('.','pt')
        fileName = path+"/psi__N_%i_L_%i_mutations_%i.csv"%(N,L,k)
    else:
        save = False

    # simTree = finiteSitesModell_investigations.simulator_KingmanFiniteSites(N,theta/2,L)
    #
    # S_redundantRowsAndColumns = simTree.getS()
    #
    # S1 = withoutNullColumns(S_redundantRowsAndColumns)
    #
    # S,n = S_and_n(S1)

    S,Nr,Nc = simData_k_mutations_total(N,L,k)

    if save:
        toCSV(S,Nr,Nc,fileName)
        print "Output saved to %s"%fileName

        # ##Test
        # Snew,nnew = psiFromCSV(fileName)
        # print "successfully loaded from %s"%fileName
        # print "S =\n",Snew,"\nn=\n",nnew

    else:
        #print "S =\n",S,"\nn=\n",n
        toCSV(S,Nr,Nc)


if __name__ == '__main__':
    main()
