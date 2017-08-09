import numpy as np
import matplotlib.pyplot as plt

from Masters.KCollisions.KNewOrbit import *

au= 149597871e3
G = 6.67408e-11
Mstar=1.2e30

a1=2.5 * au
e1= 0.998
a2= 2.4 * au
e2=0.997
'''
a1=25 * au
e1= 0.9997
a2= 25 * au
e2=0.9999
'''
I1 = 0
I2 = 10* ((2 * pi) / 360)  # in degrees

density=2000 #emily suggestioj
Qa=620
Qb=(5e-6)*density
a=0.3
b=1.5
vcol=1e5
#d2=1
def DispThres(D):
    Qd=(Qa*(D**-a))+(Qb*(D**b))
    return Qd
def DispShattering(D):
    Qs = (Qa * (D ** -a))
    return Qs
def DispGraph():
    N=100
    D=np.linspace(-1,3,N)
    Q=np.zeros((N))
    for i in range(0,N):
        D[i]=10**D[i]
    Qd=DispThres(D)
    Qs=DispShattering(D)
    Q=0.5*(vcol**2.)#0.5*((d2/D)**3.)*(vcol**2.) Even mass fragments
    plt.figure()
    plt.loglog(D,Qd,label=' Qd')
    plt.loglog(D, Qs, label=' Qs')
    #plt.loglog([D[0],D[N-1]],[Q,Q],label='Spec Energy')
    plt.xlabel('Diameter D  meters')
    plt.ylabel('Qd and Qs  J/kg')
    plt.legend()
    plt.title('Dispersion and Shattering Threshold against Diamter')
    plt.show()
    return
def flrfunc(vcol,D1,D2):
    Q = 0.5 * (vcol ** 2.)*((D2/D1)**3.)
    Qd=DispThres(D1)
    if Q<Qd: # cratering
        #print('Cratering')
        flr=1-0.5*(Q/Qd)
    if Q>=Qd: #shattering
        #print('Shattering')
        flr=0.5*((Qd/Q)**1.24)


    return flr


def ReductionGraph():
    N = 1000  # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)  # NOTE CHECK L def!! + or -
    R = np.zeros((2, N))  # radial collision points
    C = np.zeros((2, N))  # theta collision points

    for i in range(1, N):
        [[C[0, i], C[1, i]], [R[0, i], R[1, i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    Rdot1 = rdot(a1, e1, C[:, :], Mstar)
    Thetadot1 = thetadot(a1, e1, C[:, :], Mstar)
    Rdot2 = rdot(a2, e2, C[:, :] - L, Mstar)
    Thetadot2 = thetadot(a2, e2, C[:, :] - L, Mstar)
    Vinc1 = np.zeros((2, 3, N))
    Vinc2 = np.zeros((2, 3, N))
    vrelpar = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            Vinc1[x, :, i] = np.array(
                [Rdot1[x, i], R[x, i] * Thetadot1[x, i] * cos(I1), R[x, i] * Thetadot1[x, i] * sin(I1)])
            Vinc2[x, :, i] = np.array(
                [Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])
            vrelpar[x, i] = np.linalg.norm(Vinc1[x, :, i] - Vinc2[x, :, i])

    '''
    #vrel plot
    plt.figure()
    plt.semilogy(L[1:N-1],vrelpar[0,1:N-1],label='Vrel at a')
    plt.legend()
    '''

    #plt.rc('text', usetex=True)
    plt.rcParams.update({'font.size': 12})

    plt.figure()

    flr = np.zeros(N)
    dflr=np.zeros(N)
    d2=1
    for dfrac in range(-3,4):
        d1=(10**dfrac)*d2

        for i in range(1, N):
            flr[i]=flrfunc(vrelpar[0, i],d1,d2)
            dflr[i]=(flr[i])**(1/3)

        plt.semilogy(L[1:N-1]/pi,dflr[1:N-1],label='D1/D2=10^%s'%((dfrac)))


    plt.legend()
    plt.title('Reduction in D1 against '+r"$\lambda$")
    plt.xlabel(r"$\lambda$" + '/pi')
    plt.ylabel('Reduction in D1')
    plt.show()
    return

ReductionGraph()
#DispGraph()




