import numpy as np
import matplotlib.pyplot as plt
density=2700 #emily suggestioj
Qa=620
Qb=(5e-6)*density
a=0.5
b=1.5
vcol=1e5
d2=1
def DispThres(D):
    Qd=(Qa*(D**-a))+(Qb*(D**b))
    return Qd

def DispGraph():
    N=100
    D=np.linspace(-20,12,N)
    Q=np.zeros((N))
    for i in range(0,N):
        D[i]=10**D[i]
    Qd=DispThres(D)
    Q=0.5*(vcol**2.)#0.5*((d2/D)**3.)*(vcol**2.)
    #print('Q',Q)
    plt.figure()
    plt.loglog(D,Qd,label=' Dispersion')
    #plt.loglog([D[0],D[N-1]],[Q,Q],label='Spec Energy')
    plt.xlabel('Diameter D')
    plt.ylabel('Dispesion Thresh')
    plt.legend()
    plt.show()
    return

DispGraph()