
import finiteSitesModell_investigations as fsm
import time

Ls = (10,50,100,200,500)
ns = (2,10)

t0 = time.time()
for L in Ls:
    for n in ns:
        for t_factor in (0.1,1.0,10):
            steps = 20
            N = 1000
            thetaMax = t_factor*L
            print "generating plots for (n,L,thetaMax,steps,N) = (%i,%i,%4.2f,%i,%i)"%(n,L,thetaMax,steps,N)
            t1=time.time()
            fsm.generate_plot_1(n,L,thetaMax = t_factor*L,steps=steps,N = 1000)
            t2 = time.time()
            print "Time elapsed:\n\t%.3f\t(total)\n\t%.3f\t(last iteration)\n"%(t2-t0,t2-t1)
print "done!"
