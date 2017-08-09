import matplotlib.pyplot as plt
import numpy as np

from Masters.KCollisions.KNewOrbit import *

au= 149597871e3
G = 6.67408e-11
Mstar=1.2e30

a1=2.5 * au
e1= 0.998 #0.998
a2= 2.4 * au
e2=0.997 #0.995
'''
a1=25 * au
e1= 0.9997
a2= 25 * au
e2=0.9999
'''
I1 = 0
I2 = 1* ((2 * pi) / 360)  # in degrees

def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)*(10**6.) #in Yrs
    wp = (2 * np.pi) / Tp
    return (Tp,wp)

(Tp1,wp1)=TPrecess(a1,e1)
(Tp2,wp2)=TPrecess(a2,e2)

def Segments(Indexes,Function): # finds the sections for which a graph is negative, given crossing points
    #print('Indexes',Indexes)

    NRanges=np.array([[0],[0]]) #creates the starting array, delete first after
    for i in range (0,np.size(Indexes)):  #note doesnt do the last, not that it matters
        if (Function[Indexes[i]])>0:
            #print('segment added')
            NRanges=np.concatenate((NRanges,Indexes[i:i+2]) ,axis=1)
    NRanges=np.delete(NRanges,0,1)
    return NRanges

def ICLBeste(a1,e1,a2,e2,Mstar):
    N=1000 # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)   #NOTE CHECK L def!! + or -
    R= np.zeros((2,N))  #radial collision points
    C = np.zeros((2, N))   #theta collision points

    for i in range(0,N):
        [[C[0,i], C[1,i]], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    ThetaKepler=np.sqrt((G*Mstar)/(np.power((R[:,:]),3)))
    Rdot1=rdot(a1,e1,C[:,:],Mstar) #another (2,N) array. Note taken out masses, no differnece!
    Thetadot1=thetadot(a1,e1,C[:,:],Mstar)
    ThetadotK1=Thetadot1-ThetaKepler
    Rdot2 = rdot(a2, e2, C[:, :]-L,  Mstar)
    Thetadot2=thetadot(a2, e2, C[:, :]-L,  Mstar)
    ThetadotK2 = Thetadot2-ThetaKepler
    BestRdot=np.zeros((2,N))
    BestThetadot=np.zeros((2,N))
    for j in range(0,2): #two collision point loop
        for i in range(0,N): #running down all the L values
            #Rdot Values
            if np.sign(Rdot1[j,i]*Rdot2[j,i]) == 1:
                if np.array([abs(Rdot1[j,i]),abs(Rdot2[j,i])]).argmin() == 0:
                    BestRdot[j, i] = Rdot1[j,i]
                else:
                    BestRdot[j, i] = Rdot2[j, i]
            else:
                BestRdot[j,i]=0
            #ThetaDot Valeus ,
            if np.sign(ThetadotK1[j,i]*ThetadotK2[j,i]) == 1: #' closest to kepler'
                if np.array([abs(ThetadotK1[j,i]),abs(ThetadotK2[j,i])]).argmin() == 0:
                    BestThetadot[j, i] = Thetadot1[j,i]
                else:
                    BestThetadot[j, i] = Thetadot2[j, i]
            else:
                BestThetadot[j,i]=ThetaKepler[j,i]   #or exactly kepler






    Beste=newe(BestRdot, BestThetadot, R, Mstar)


    return Beste

def ICLBestePlot(a1, e1, a2, e2, Mstar): #note masses make no difference
    Beste=ICLBeste(a1, e1, a2, e2, Mstar)
    N = 1000  # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)


    plt.figure()
    plt.plot(L / pi, Beste[0, :], label='BestEccentricity of Collision a')
    plt.plot(L / pi, Beste[1, :], label='BestEccentricity of Collision b')
    plt.legend()
    plt.xlabel('Lambda/pi')
    plt.ylabel('Best Eccentricity')

    plt.show()
    return


#InitialComponents(2*au,0.99,2.1*au,0.993,1.2e30,2e10,2e10,10000)

#ICLBestePlot(2*au,0.99,2.1*au,0.993,1.2e30)


def ICLVelocitiesGraph():
    N=1000 # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)   #NOTE CHECK L def!! + or -
    R= np.zeros((2,N))  #radial collision points
    C = np.zeros((2, N))   #theta collision points

    for i in range(1,N):
        [[C[0,i], C[1,i]], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    #ThetaKepler=np.sqrt((G*Mstar)/(np.power((R[:,:]),3)))
    Rdot1=rdot(a1,e1,C[:,:],Mstar) #another (2,N) array. Note taken out masses, no differnece!
    Thetadot1=thetadot(a1,e1,C[:,:],Mstar)
    RThetadot1=R*Thetadot1
   # ThetadotK1=Thetadot1-ThetaKepler
    Rdot2 = rdot(a2, e2, C[:, :]-L,  Mstar)
    Thetadot2=thetadot(a2, e2, C[:, :]-L,  Mstar)
    RThetadot2 = R * Thetadot2
   # ThetadotK2 = Thetadot2-ThetaKepler
    '''
    Vinc1=np.zeros((2,3,N))
    Vinc2=np.zeros((2,3,N))
    vrelpar = np.zeros((2, N))
    for i in range(1,N-1):
        for x in (0,1):
            Vinc1[x,:,i]=np.array([Rdot1[x,i],R[x,i]*Thetadot1[x,i]*cos(I1),R[x,i]*Thetadot1[x,i]*sin(I1)])
            Vinc2[x, :, i] = np.array(
                [Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])

            vrelpar[x, i] = np.linalg.norm(Vinc1[x, :, i] - Vinc2[x, :, i])
    '''

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
   
    
    # graph of relative parent velocities - including Inclination
    
    # print(L[0]/pi,L[N-1]/pi)
    # print('vrelpar at L=0',vrelpar[0])
    # print('vrelpar at L=2pi', vrelpar[N-1])
    plt.rc('text', usetex=True)
    plt.rcParams.update({'font.size': 15})
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, vrelpar[0, 1:N - 2], label=r"$v_{\textnormal{col}}$" +' at A')
    plt.semilogy(L[1:N - 2] / pi, vrelpar[1, 1:N - 2], label=r"$v_{\textnormal{col}}$" +'at B')
    plt.xlabel(r"$\lambda$" +'/pi' )
    plt.ylabel(r"$v_{\textnormal{col}}$" +'m/s' )
    plt.legend()
    plt.title(r"$v_{\textnormal{col}}$")#Inclination %s degrees' % (abs(I1 - I2) * (360 / (2 * pi))))

    
    
    #h/dr
    hoverdr = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0,1):
            #a=np.array([Rdot1[x,i],R[x,i]*Thetadot1[x,i]*cos(I1),R[x,i]*Thetadot1[x,i]*sin(I1)])
            #b = np.array([Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])
            hoverdr[x, i] =(((R[x,i]**2.)*Thetadot1[x,i]*Thetadot2[x,i]*sin(I1-I2))/np.linalg.norm(np.cross(Vinc1[x,:,i],Vinc2[x,:,i])))
    hoverdr=abs(hoverdr)

    dhdt = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            # a=np.array([Rdot1[x,i],R[x,i]*Thetadot1[x,i]*cos(I1),R[x,i]*Thetadot1[x,i]*sin(I1)])
            # b = np.array([Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])
            dhdt[x, i] = (((R[x, i] ** 2.) * Thetadot1[x, i] * Thetadot2[x, i] * sin(I1 - I2)) / np.linalg.norm(
                np.cross(Vinc1[x, :, i], Vinc2[x, :, i])))*(((wp1)*(Rdot1[x,i]/Thetadot1[x,i]))-(wp2)*(Rdot2[x,i]/Thetadot2[x,i]))
    dhdt = abs(dhdt)
    '''
    #contact time
    '''
    OLD METHOD
    TcontactoRb = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            TcontactoRb[x,i]=4/(hoverdr[x,i]*abs(((wp1)*(Rdot1[x,i]/Thetadot1[x,i]))-(wp2)*(Rdot2[x,i]/Thetadot2[x,i])))
    '''
    '''
    TcontactoRb = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            TcontactoRb[x, i] = 4 / dhdt[x,i]





    #hoverdrGraph
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, hoverdr[0, 1:N - 2], label='h/dr at a')
    plt.semilogy(L[1:N - 2] / pi, hoverdr[1, 1:N - 2], label='h/dr at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('dh/dr, dimensionless')
    plt.legend()

    #dhdtGraph
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, dhdt[0, 1:N - 2], label='dh/dt at a')
    plt.semilogy(L[1:N - 2] / pi, dhdt[1, 1:N - 2], label='dh/dt at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('dh/dt, m/yr')
    plt.legend()

    #ContactTime
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, TcontactoRb[0, 1:N - 2], label='Tcontact/Rbeam at a')
    plt.semilogy(L[1:N - 2] / pi, TcontactoRb[1, 1:N - 2], label='Tcontact/Rbeam at b')
    plt.xlabel('Lambda/pi')
    plt.ylabel('Tcontact/Rbeam, years')
    plt.title('Inclination difference of %s degrees'%(abs(I1-I2)*(360/(2*pi))))
    plt.legend()

    '''
    #low inclination
    LowTcontactoRb= np.zeros((2, N))
    # *(I2-I1) Removed I dep
    for i in range(1, N - 1):
        for x in (0, 1):
            LowTcontactoRb[x,i]=(4*((Rdot1[x,i]*Thetadot2[x,i])-(Rdot2[x,i]*Thetadot1[x,i])))/((((wp1)*(Rdot1[x,i]/Thetadot1[x,i]))-(wp2)*(Rdot2[x,i]/Thetadot2[x,i]))*R[x,i]*Thetadot1[x,i]*Thetadot2[x,i])
    LowTcontactoRb=abs(LowTcontactoRb)
    #low inclination
    dhdt = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            # a=np.array([Rdot1[x,i],R[x,i]*Thetadot1[x,i]*cos(I1),R[x,i]*Thetadot1[x,i]*sin(I1)])
            # b = np.array([Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])
            dhdt[x, i] = 4/LowTcontactoRb[x,i]
    dhdt = abs(dhdt)


    '''
    plt.rc('text', usetex=True)
    plt.rcParams.update({'font.size': 12})
    # dhdtGraph
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, dhdt[0, 1:N - 2],  label = r"$\frac{1}{\delta I}\frac{\delta h}{\delta t}$"+' at A')
    plt.semilogy(L[1:N - 2] / pi, dhdt[1, 1:N - 2],  label = r"$\frac{1}{\delta I}\frac{\delta h}{\delta t}$"+' at B')
    plt.xlabel(r"$\lambda$" + '/pi')
    plt.ylabel(r"$\frac{1}{\delta I}\frac{\delta h}{\delta t}$"+' m/year')
    plt.legend()
    '''
    plt.rc('text', usetex=True)
    plt.rcParams.update({'font.size': 14})
    # Tcontact/IRb
    '''

    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, LowTcontactoRb[0, 1:N - 2], label= 'value at A')
    plt.semilogy(L[1:N - 2] / pi, LowTcontactoRb[1, 1:N - 2], label='value at B')
    plt.xlabel(r"$\lambda$" + '/pi')
    plt.ylabel(r"$\frac{1}{R_{\small{Beam}}\delta I}T_{\small{Contact}}$" + ' year/m')
    plt.title(r"$\frac{1}{R_{\small{Beam}}\delta I}T_{\small{Contact}}$" + ' against '+r"$\lambda$" + '/pi' )
    plt.legend()
    '''
    LowTcontactoRb=(10**5.)*LowTcontactoRb
    plt.figure()
    plt.semilogy(L[1:N - 2] / pi, LowTcontactoRb[0, 1:N - 2], label= 'value at A')
    plt.semilogy(L[1:N - 2] / pi, LowTcontactoRb[1, 1:N - 2], label='value at B')
    plt.xlabel(r"$\lambda$" + '/pi')
    plt.ylabel(r"$T_{\small{Contact}}$" + ' year')
    plt.title(r"$T_{\small{Contact}}$" + ' against '+r"$\lambda$" + '/pi' )
    plt.legend()

    plt.show()


    return



ICLVelocitiesGraph()