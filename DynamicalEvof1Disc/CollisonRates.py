import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from KCollisions.KNewOrbit import *

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
I2 = 10 * ((2 * pi) / 360)  # in degrees


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
    for x in range(0, Xnumber):
        xdepIntdata[x]=(IntTimeData[x]**2.)/IntVolData[x]
    return xdepIntdata


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
    plt.show()
    return





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


    T1=2*pi*np.sqrt((a1**3.)/(G*Mstar))
    T2 = 2 * pi * np.sqrt((a2 ** 3.) / (G * Mstar))

    CollisionRateData=np.zeros((2,N))
    for i in range(1, N - 1):
        for x in (0, 1):
            CollisionRateData[x,i]=(1/(2*pi))*(1/(T1*T2))*(vrelpar[x,i]/(v1[x,i]*v2[x,i]*sin(gammangle[x,i])))#

    '''
    plt.figure()
    plt.plot(L[1:N - 1] / pi, 1/np.sin(gammangle[0,1:N-1]), label='1/sin(gamma) at a')
    plt.plot(L[1:N - 1] / pi, 1/np.sin(gammangle[1,1:N-1]), label='1/sin(gamma) at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('1/sin(gamma)')
    plt.legend()
    '''
    '''
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

    plt.figure()
    plt.semilogy(L[1:N-1]/pi,CollisionRateData[0,1:N-1],label='g at a')
    plt.semilogy(L[1:N - 1]/pi, CollisionRateData[1, 1:N - 1], label='g at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('g(lambda)')
    plt.title('g(lambda), dependance of the Collision Rate on orbits')
    plt.legend()
    plt.show()
    return

#CollisionRateCalc()




IntegralGraph()

'''
 #graph of all velocities
 plt.figure()
 plt.subplot(211)
 plt.plot(L/(pi),Rdot1[0,:], label='rd1 at a')
 plt.plot(L/(pi),RThetadot1[0,:],label='rtd1 at a')
 plt.plot(L/(pi), Rdot2[0, :], label='rd2 at a')
 plt.plot(L / ( pi), RThetadot2[0, :], label='rtd2 at a')
 plt.legend()
 plt.subplot(212)
 plt.plot(L / ( pi), Rdot1[1, :], label='rd1 at b')
 plt.plot(L / ( pi), RThetadot1[1, :], label='rtd1 at b')
 plt.plot(L / ( pi), Rdot2[1, :], label='rd2 at b')
 plt.plot(L / (pi), RThetadot2[1, :], label='rtd2 at b')
 plt.legend()

 #graph of angles
 plt.figure()
 angle1=np.arctan(Rdot1/RThetadot1)
 angle2 = np.arctan(Rdot2 / RThetadot2)
 plt.subplot(211)
 plt.plot(L/pi,abs(angle1[0,:]/pi),label='angle 1 at a')
 plt.plot(L / pi, abs(angle2[0, :]/pi), label='angle 2 at a')
 plt.legend()
 plt.subplot(212)
 plt.plot(L / pi, abs(angle1[1, :]/pi), label='angle 1 at b')
 plt.plot(L / pi, abs(angle2[1, :]/pi), label='angle 2 at b')
 plt.legend()

 #graph of relative parent velocities
 vrelpar=np.zeros((2,N))
 for i in range(1,N-1):
     vrelpar[:,i]=np.sqrt(((Rdot1[:,i]-Rdot2[:,i])**2.)+((R[:,i]*(Thetadot1[:,i]-Thetadot2[:,i]))**2.))
 #print(L[0]/pi,L[N-1]/pi)
 #print('vrelpar at L=0',vrelpar[0])
 #print('vrelpar at L=2pi', vrelpar[N-1])
 plt.figure()
 plt.semilogy(L[1:N-2]/pi,vrelpar[0,1:N-2],label='vrelpar at a')
 plt.semilogy(L[1:N - 2] / pi, vrelpar[1, 1:N - 2], label='vrelpar at b')
 plt.xlabel('L/pi')
 plt.ylabel('Relative Velocities of the parents, m/s')
 plt.legend()
 '''
