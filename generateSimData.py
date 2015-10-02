from almostInfiniteSites_simulateData import simData,toCSV
from time import ctime

path = "./simData"

n_list = [2**i for i in range(1,11)]
#theta_list = [10.0**i for i in range(-2,4)]
theta_list = [10.0**2]
L_list = [1,2,4,8,16,32,64,128,256,512,1024,2048]

arglist = [(a,b,c) for a in n_list for b in theta_list for c in L_list]

for n,theta,L in arglist:
    thetaStr = ("%1.4f"%theta).replace('.','pt')
    fileName = path+"/psi__N_%i_theta_%s_L_%i.csv"%(n,thetaStr,L)

    print "%s\tsimulating data for n,theta,L = %i,%.2f,%i"%(ctime(),n,theta,L)
    S,n = simData(n,theta,L)
    toCSV(S,n,fileName)
    print "saved to %s\n"%fileName
