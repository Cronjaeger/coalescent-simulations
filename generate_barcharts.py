# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:31:49 2015

@author: mathias
"""

import finiteSitesModell_investigations

if __name__ == "__main__":
    arglist = []
    for L in [2, 8, 32, 128, 512, 2048, 8192]:
#    for L in [2,8,32]:
        for n in (2,8,32,128):
            arglist.append([1000,L,n])

#    print arglist

    finiteSitesModell_investigations.run_generatePlot_of_mutationTypes(arglist)