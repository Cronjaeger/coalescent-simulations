# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 08:58:43 2015

@author: mathias
"""

import finiteSitesModell_investigations

if __name__ == "__main__":
    arglist = []
    for L in [2, 8, 32, 128, 512, 2048, 8192]:
        for n in (2,8,32,128):
            arglist.append([1000,L,n])

    finiteSitesModell_investigations.run_generateScatterplots(arglist)