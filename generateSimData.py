from almostInfiniteSites_simulateData import simData,toCSV,simData_k_mutations_total
from time import ctime


# simulation_type = 'regular'
simulation_type = 'fixed number of segregating sites'
EXCESS_MUTATIONS_MUST_BE_VISIBLE = True

path = "./simData/fixed_number_of_segregating_sites"

# n_list = [2**i for i in range(1,11)]
# theta_list = [10.0**i for i in range(-2,4)]
# L_list = [1,2,4,8,16,32,64,128,256,512,1024,2048]
# n_excess_mutations_list = [1,2,3,4,5]

n_list = [i for i in range(2,20)]
k_list = [[i for i in range(1,20)]]
L = 1000
theta = 100

## This is outdated, and relies on old .csv-encoding
# if simulation_type == 'regular':
#
#     arglist = [(a,b,c) for a in n_list for b in theta_list for c in L_list]
#     path = "./simData"
#
#     for n,theta,L in arglist:
#         thetaStr = ("%1.4f"%theta).replace('.','pt')
#         fileName = path+"/psi__N_%i_theta_%s_L_%i.csv"%(n,thetaStr,L)
#
#         print "%s\tsimulating data for n,theta,L = %i,%.2f,%i"%(ctime(),n,theta,L)
#         S,n = simData(n,theta,L)
#         toCSV(S,n,fileName)
#         print "saved to %s\n"%fileName

if simulation_type == 'fixed number of segregating sites':
    n_list = [i for i in range(2,21)]
    k_list = [i for i in range(0,21)]
    L = 1000
    samples = range(1,31)
    # theta = 100
    argList = [(a,b,c) for a in n_list for b in k_list for c in samples]

    for n,k,sample_no in argList:
        # print n,k,L
        fileName = path+"/psi__N_%i_L_%i__mutations_%i__sample_%i.csv"%(n,L,k,sample_no)
        S,Nr,Nc = simData_k_mutations_total(n,L,k)

        # while len(Nc) < k:
        #     S,Nr,Nc = simData_k_mutations_total(n,L,k)

        toCSV(S,Nr,Nc,fileName)
        print "saved to %s\n"%fileName


else:
    print "unknown type of simulation '%s', sorry..."%simulation_type
    # argList = [(n,k,L) for k in n_excess_mutations_list for n in n_list for L in L_list]
    #
    # if EXCESS_MUTATIONS_MUST_BE_VISIBLE:
    #     path = "./simData/conditioned/visibleExces"
    #     fileNameDescription = "visibleExces"
    # else:
    #     path = "./simData/conditioned/fixedExces"
    #     fileNameDescription = "fixedExces"
