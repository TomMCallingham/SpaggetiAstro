import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from Masters.KCollisions.KNewOrbit import *

au = 149597871e3
G = 6.67408e-11
Mstar = 1.2e30

a1 = 2.5 * au
e1 = 0.998
a2 = 2.4 * au
e2 = 0.997
'''
a1=25 * au
e1= 0.9997
a2= 25 * au
e2=0.9999
'''
I1 = 0
I2 = 1* ((2 * pi) / 360)  # in degrees


def TPrecess(a, e):
    Tp = 0.15 * ((1 - (e ** 2.)) / (1 - (0.999 ** 2.))) * ((a / au) ** 2.5) * (10 ** 6.)  # in Yrs
    wp = (2 * np.pi) / Tp
    return (Tp, wp)


(Tp1, wp1) = TPrecess(a1, e1)
(Tp2, wp2) = TPrecess(a2, e2)


def IntegrandVol(z,x):
    I=np.sqrt((1-(z**2.))*(1-((z-x)**2.)))
    return I
def IntVol(Xnumber):
    X=np.linspace(0,2,Xnumber)
    IntVolData=np.zeros(Xnumber)
    for x in range(0,Xnumber):
        IntVolData[x]=quad(IntegrandVol,X[x]/2,1, args=(X[x]))[0]

    return IntVolData
def IntTime(Xnumber):
    X = np.linspace(0, 2, Xnumber)
    IntTimeData = np.zeros(Xnumber)
    for x in range(0, Xnumber):
        IntTimeData[x] = quad(IntegrandVol, X[x]-1 , 1, args=(X[x]))[0]

    return IntTimeData

def xdepInt(Xnumber):
    X = np.linspace(0, 2, Xnumber)
    IntTimeData = IntTime(Xnumber)
    IntVolData=IntVol(Xnumber)
    xdepIntdata=np.zeros(Xnumber)
    for x in range(0, Xnumber-1):
        xdepIntdata[x]=(3/8)*(IntTimeData[x]**2.)/IntVolData[x]
    return xdepIntdata




def IntAverargeXdep():
    Xnumber=1000
    X = np.linspace(0, 2, Xnumber)
    xdepIntData = xdepInt(Xnumber)
    AverargeXdep = 0
    for x in range(0,Xnumber-1):
        AverargeXdep=AverargeXdep+xdepIntData[x]
    AverargeXdep=(2/Xnumber)*AverargeXdep
    print('AverageXdep',AverargeXdep)
    return


def IntegralGraph():
    Xnumber=1000
    IntVolData=IntVol(Xnumber)
    xdepIntData=xdepInt(Xnumber)
    X = np.linspace(0, 2, Xnumber)
    plt.figure()
    plt.plot(X,IntVolData)
    plt.xlabel('x=h/Rbeam')
    plt.ylabel('IntVol')
    plt.figure()
    plt.plot(X, xdepIntData)
    plt.xlabel('x=h/Rbeam')
    plt.ylabel('f(x)')
    plt.title('Rcol dependance on x=h/Rbeam')
    print('finalised graphing')

    return


#IntAverargeXdep()
#IntegralGraph()




def CollisionRateCalc():
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
    v1=np.zeros((2, N))
    v2=np.zeros((2, N))
    gammangle=np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            Vinc1[x, :, i] = np.array(
                [Rdot1[x, i], R[x, i] * Thetadot1[x, i] * cos(I1), R[x, i] * Thetadot1[x, i] * sin(I1)])
            Vinc2[x, :, i] = np.array(
                [Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])
            vrelpar[x, i] = np.linalg.norm(Vinc1[x, :, i] - Vinc2[x, :, i])
            v1[x,i]=np.sqrt((Rdot1[x,i]**2.)+((R[x,i]*Thetadot1[x,i])**2.))
            v2[x, i] = np.sqrt((Rdot2[x, i] ** 2.) + ((R[x, i] * Thetadot2[x, i] )** 2.))
            gammangle[x,i]=np.arccos(np.vdot(Vinc1[x,:,i],Vinc2[x,:,i])/(v1[x,i]*v2[x,i]))


    T1=2*pi*np.sqrt((a1**3.)/(G*Mstar)) #Orbital Period
    T2 = 2 * pi * np.sqrt((a2 ** 3.) / (G * Mstar))

    gCollisionRateData=np.zeros((2,N))
    for i in range(1, N - 1):
        for x in (0, 1):
            gCollisionRateData[x,i]=(16/(3*pi))*(1/(T1*T2))*(vrelpar[x,i]/(v1[x,i]*v2[x,i]*sin(gammangle[x,i])))#

    '''
    plt.figure()
    plt.plot(L[1:N - 1] / pi, 1/np.sin(gammangle[0,1:N-1]), label='1/sin(gamma) at a')
    plt.plot(L[1:N - 1] / pi, 1/np.sin(gammangle[1,1:N-1]), label='1/sin(gamma) at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('1/sin(gamma)')
    plt.legend()

    plt.figure()
    plt.plot(L[1:N-1]/pi,v1[0,1:N-1],label='v1 at a')
    plt.plot(L[1:N - 1]/pi, v2[0, 1:N - 1], label='v2 at a')
    plt.xlabel('Lambda/pi')
    plt.ylabel('v')
    plt.legend()

    plt.figure()
    plt.plot(L[1:N-1]/pi,gammangle[0,1:N-1],label='gammange at a')
    plt.plot(L[1:N - 1]/pi, gammangle[1,1:N-1], label='gammangle at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('v')
    plt.legend()
    '''
    plt.rc('text', usetex=True)
    plt.rcParams.update({'font.size': 15})

    plt.figure()
    plt.semilogy(L[1:N-1]/pi,gCollisionRateData[0,1:N-1]*year,label='g at A')
    plt.semilogy(L[1:N - 1]/pi, gCollisionRateData[1, 1:N - 1]*year, label='g at B')
    plt.xlabel(r"$\lambda$" + '/pi')
    plt.ylabel('g /years /meters ')
    plt.title('function g against ' + r"$\lambda$") #dependance of the Collision Rate on orbits against lambda. I=%s%((I2-I1)*(360/(2*pi))))
    plt.legend()

    return

CollisionRateCalc()



#IntAverargeXdep()
#IntegralGraph()

def AreaDep():
    D1=np.linspace(1,100,1000)
    ncross=np.zeros((1000))
    D2=[1,5,10,50,100]
    plt.figure()
    K=2*(10**11.)
    plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    plt.rcParams.update({'font.size':12})
    for d2 in range(0,5):

        for d1 in range(0,1000):
            ncross[d1]=(pi*((D1[d1]+D2[d2])**2.)*(D1[d1]**(-3.5))/4)
        #plt.semilogy(D1,K*ncross,label=r"$ D_{2}=$"+'%s m'%D2[d2])
        #
    plt.legend()
    plt.xlabel(r"$D_{1}$")
    plt.ylabel(r"$n\left(D_{1}\right)\sigma_{12}$"+' / '+r"$\textnormal{m}^{2}$")
    plt.title(r"$n\left(D_{1}\right)\sigma_{12}$"+' Against '+r"$D_{1}$"+' for different ' +r"$D_{2}$")



    return

plt.show()
#AreaDep()