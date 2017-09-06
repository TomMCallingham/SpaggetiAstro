from scipy.special import hyp2f1
import numpy as np
from Functions import *


c=3e8
density=2e3
Lsun=3.828e26
Tini=0.5e9
L=3.26*Lsun*((0.1+Tini*(10**-6.))**(-1.18))
Q=1
alpha=((3*L*Q)/(16*np.pi*(c**2.)*density))
beta=alpha/(G*Mstar)
#print('Beta',beta)

def TPrAnalytic(a,e):

    TPr=(((2*((a*(1-(e**2.)))**2.))/(5*alpha))*((1/np.sqrt(1-(e**2.)))-((3/8)*hyp2f1(0.5, 0.8, 1.8,e**2))))/year
    return TPr




def TprGraphsae():
    print('creating Tpr data...')
    N = 1000
    e=np.linspace(0,0.99,N)
    a=np.linspace(0,10,N)*au
    Tpr=np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            Tpr[i,j]=TPrAnalytic(a[i],e[j])
    print('plotting...')
    plt.figure()
    evals=[0.1,0.5,0.8,0.9,0.99]
    for i in range(0,np.size(evals)):
        evalsn=int(evals[i]/(e[1]-e[0]))
        plt.loglog(a/au,Tpr[:,evalsn], label='e=%.2g'%evals[i])

    plt.xlabel('Semi Major Axis au')
    plt.ylabel('Tpr, yrs')
    plt.title('Tpr vs a')
    plt.legend()

    plt.figure()
    print('plotted')
    avals=np.arange(1,16,2)*au
    for i in range(1, np.size(avals)):
        avalsn=int(avals[i]/(avals[1]-avals[0]))
        plt.semilogy(e, Tpr[avalsn,:], label='a=%.2g au' % (avals[i]/au))

    plt.xlabel('e')
    plt.ylabel('Tpr, yrs')
    plt.title('Tpr vs e')
    plt.legend()

    return


#TprGraphsae()

plt.show()
#print(np.log10(TPrAnalytic(10*au,0.9993)))
def TPrhighelim():
    q=0.01 #in au
    epsilon=(32*pi*(c**2.)*density)/(15*L*Q)
    constant=2*epsilon*(au**2.)*(q**(3/2))/year
    print(constant)
    return

#TPrhighelim()
