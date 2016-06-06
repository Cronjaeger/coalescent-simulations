# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 17:31:49 2015

@author: mathias
"""

import finiteSitesModell_investigations

N_events_to_Stop_after = 1

if __name__ == "__main__":
    arglist = []
    #for L in [2, 8, 32, 128, 512, 2048, 8192]:
    for L in [8,64,128,512,2048]:
        for n in [8,32,128]:
            arglist.append([1000,L,n])

    # Add some BIG parameters (in case the computer terminated too early)
    arglist.append([1000,2048,8])
    arglist.append([1000,2048,128])
#    print arglist

    finiteSitesModell_investigations.run_generatePlot_of_mutationTypes(arglist,N_events_to_Stop_after)
