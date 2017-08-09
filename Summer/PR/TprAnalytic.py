from scipy.special import hyp2f1
import numpy as np
from Functions import *
Lsun=3.828e26
#Values
a1 = 2.5 * au
e1 = 0.998
a2 = 2.4 * au
e2 = 0.997
I=1* ((2 * pi) / 360)
c=3e8
density=2e3
Lsun=3.828e26
Tini=0.5e9
L=3.26*Lsun*((0.1+Tini*(10**-6.))**(-1.18))
Q=1
alpha=((3*L*Q)/(16*np.pi*(c**2.)*density))

def TPrAnalytic(a,e):

    TPr=(((2*((a*(1-(e**2.)))**2.))/(5*alpha))*((1/np.sqrt(1-(e**2.)))-((3/8)*hyp2f1(0.5, 0.8, 1.8,e**2))))/year
    return TPr




def TprGraphs():
    print('creating data...')
    N = 1000
    e=np.linspace(0,0.9,N)
    a=np.linspace(0,5,N)*au
    Tpr=np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            Tpr[i,j]=TPrAnalytic(a[i],e[j])
    print('plotting...')
    plt.figure()
    J = 5
    for i in range(1,J):
        plt.loglog(a/au,Tpr[:,int(N*(i/J))-1], label='e=%.2g'%e[int(N*(i/J))-1])

    plt.xlabel('Semi Major Axis au')
    plt.ylabel('Tpr, yrs')
    plt.title('Tpr vs a')
    plt.legend()

    plt.figure()
    print('plotted')

    for i in range(1, J):
        plt.semilogy(e, Tpr[int(N*(i/J))-1,:], label='a=%.2g au' % (a[int(N*(i/J))-1]/au))

    plt.xlabel('e')
    plt.ylabel('Tpr, yrs')
    plt.title('Tpr vs e')
    plt.legend()

    return

#TprGraphs()

plt.show()
#print(np.log10(TPrAnalytic(10*au,0.9993)))