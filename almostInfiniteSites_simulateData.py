# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 08:36:25 2015

@author: mathias
"""

import finiteSitesModell_investigations
import sys
import numpy as np

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

def toCSV(S,n,fileName):
    a = np.c_[S,n]

    #sort rows so that n is in non-ascending order
    a = a[a[:,-1].argsort()[::-1]]

    np.savetxt(fileName, a, fmt = '%d', delimiter=", ")

def main():
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

    simTree = finiteSitesModell_investigations.simulator_KingmanFiniteSites(N,theta/2,L)

    S_redundantRowsAndColumns = simTree.getS()

    S1 = withoutNullColumns(S_redundantRowsAndColumns)

    S,n = S_and_n(S1)

    if save:
        toCSV(S,n,fileName)
        print "Output saved to %s"%fileName

    else:
        print "S =\n",S,"\nn=\n",n

if __name__ == '__main__':
    main()
