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
e2 = 0.995
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



def ICLVelocitiesGraph():
    N = 1000  # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)  # NOTE CHECK L def!! + or -
    R = np.zeros((2, N))  # radial collision points
    C = np.zeros((2, N))  # theta collision points

    for i in range(1, N):
        [[C[0, i], C[1, i]], [R[0, i], R[1, i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    # ThetaKepler=np.sqrt((G*Mstar)/(np.power((R[:,:]),3)))
    Rdot1 = rdot(a1, e1, C[:, :], Mstar)  # another (2,N) array. Note taken out masses, no differnece!
    Thetadot1 = thetadot(a1, e1, C[:, :], Mstar)
    RThetadot1 = R * Thetadot1
    # ThetadotK1=Thetadot1-ThetaKepler
    Rdot2 = rdot(a2, e2, C[:, :] - L, Mstar)
    Thetadot2 = thetadot(a2, e2, C[:, :] - L, Mstar)
    RThetadot2 = R * Thetadot2
    # ThetadotK2 = Thetadot2-ThetaKepler

    vinc1 = np.zeros((2, 3, N))
    vinc2 = np.zeros((2, 3, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            vinc1[x, :, i] = np.array(
                [Rdot1[x, i], R[x, i] * Thetadot1[x, i] * cos(I1), R[x, i] * Thetadot1[x, i] * sin(I1)])
            vinc2[x, :, i] = np.array(
                [Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])

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

    # graph of relative parent velocities - including Inclination
    vrelpar = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            vrelpar[x, i] = np.linalg.norm(vinc1[x, :, i] - vinc2[x, :, i])
    # print(L[0]/pi,L[N-1]/pi)
    # print('vrelpar at L=0',vrelpar[0])
    # print('vrelpar at L=2pi', vrelpar[N-1])
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, vrelpar[0, 1:N - 2], label='vrelpar at a')
    plt.semilogy(L[1:N - 2] / pi, vrelpar[1, 1:N - 2], label='vrelpar at b')
    plt.xlabel('L/pi')
    plt.ylabel('Relative Velocities of the parents, m/s')
    plt.legend()
    plt.title('Relative Velocities of the parents, Inclination %s degrees' % (abs(I1 - I2) * (360 / (2 * pi))))

    # h/dr
    hoverdr = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            # a=np.array([Rdot1[x,i],R[x,i]*Thetadot1[x,i]*cos(I1),R[x,i]*Thetadot1[x,i]*sin(I1)])
            # b = np.array([Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])
            hoverdr[x, i] = (((R[x, i] ** 2.) * Thetadot1[x, i] * Thetadot2[x, i] * sin(I1 - I2)) / np.linalg.norm(
                np.cross(vinc1[x, :, i], vinc2[x, :, i])))
    hoverdr = abs(hoverdr)

    # contact time
    TcontactoRb = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            TcontactoRb[x, i] = 4 / (
            hoverdr[x, i] * abs(((wp1) * (Rdot1[x, i] / Thetadot1[x, i])) - (wp2) * (Rdot2[x, i] / Thetadot2[x, i])))

    AltTcontactoRb = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            AltTcontactoRb[x, i] = (4 * ((Rdot1[x, i] * Thetadot2[x, i]) - (Rdot2[x, i] * Thetadot1[x, i]))) / (
            (((wp1) * (Rdot1[x, i] / Thetadot1[x, i])) - (wp2) * (Rdot2[x, i] / Thetadot2[x, i])) * R[x, i] * Thetadot1[
                x, i] * Thetadot2[x, i] * (I2 - I1))
    AltTcontactoRb = abs(AltTcontactoRb)

    # hoverdrGraph
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, hoverdr[0, 1:N - 2], label='h/dr at a')
    plt.semilogy(L[1:N - 2] / pi, hoverdr[1, 1:N - 2], label='h/dr at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('h/dr, dimensionless')
    plt.legend()

    # ContactTime
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, TcontactoRb[0, 1:N - 2], label='Tconact/Rbeam at a')
    plt.semilogy(L[1:N - 2] / pi, TcontactoRb[1, 1:N - 2], label='Tconact/Rbeam at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('Tconact/Rbeam, years')
    plt.title('Inclination difference of %s degrees' % (abs(I1 - I2) * (360 / (2 * pi))))
    plt.legend()

    # AltContactTime
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, AltTcontactoRb[0, 1:N - 2], label='Tconact/Rbeam at a')
    plt.semilogy(L[1:N - 2] / pi, AltTcontactoRb[1, 1:N - 2], label='Tconact/Rbeam at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('AltTconact/Rbeam, years')
    plt.title('Inclination difference of %s degrees' % (abs(I1 - I2) * (360 / (2 * pi))))
    plt.legend()

    plt.show()

    return


'''
def IntegrandCross(z):
    I=np.sqrt(1-(z**2.))
    return I
def IntCross(Xnumber):
    X=np.linspace(0,2,Xnumber)
    IntCrossData=np.zeros(Xnumber)
    for x in range(0,Xnumber):
        IntCrossData[x]=quad(IntegrandCross,X[x]-1,1)[0]
    return IntCrossData
'''

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
'''
#v1
def xdepInt(Xnumber):
    X = np.linspace(0, 2, Xnumber)
    IntCrossData = IntCross(Xnumber)
    IntVolData=IntVol(Xnumber)
    xdepIntdata=np.zeros(Xnumber)
    for x in range(0, Xnumber):
        xdepIntdata[x]=(IntCrossData[x]**4.)/(((2-X[x])**2.)*IntVolData[x])
    return xdepIntdata
#v2
def xdepInt(Xnumber):
    X = np.linspace(0, 2, Xnumber)
    IntCrossData = IntCross(Xnumber)
    IntVolData=IntVol(Xnumber)
    xdepIntdata=np.zeros(Xnumber)
    for x in range(0, Xnumber):
        xdepIntdata[x]=(IntCrossData[x]**2.)/(IntVolData[x])
    return xdepIntdata
'''
def xdepInt(Xnumber):
    X = np.linspace(0, 2, Xnumber)
    IntTimeData = IntTime(Xnumber)
    IntVolData=IntVol(Xnumber)
    xdepIntdata=np.zeros(Xnumber)
    for x in range(0, Xnumber):
        xdepIntdata[x]=(IntTimeData[x]**2.)/IntVolData[x]
    return xdepIntdata


def IntegralGraph(Xnumber):
    IntVolData=IntVol(Xnumber)
    '''
    #IntCrossData=IntCross(Xnumber)
    plt.figure()
    plt.plot(X, IntCrossData)
    plt.xlabel('x=h/Rbeam')
    plt.ylabel('IntCross')
    '''
    xdepIntData=xdepInt(Xnumber)
    X = np.linspace(0, 2, Xnumber)
    plt.figure()
    plt.plot(X,IntVolData)
    plt.xlabel('x=h/Rbeam')
    plt.ylabel('IntVol')



    plt.figure()
    plt.plot(X, xdepIntData)
    plt.xlabel('x=h/Rbeam')
    plt.ylabel('xdepInt')




    print('finised graphinh')
    plt.show()

    return

IntegralGraph(100)

