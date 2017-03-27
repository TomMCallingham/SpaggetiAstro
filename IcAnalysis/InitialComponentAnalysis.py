import matplotlib.pyplot as plt
import numpy as np

from KCollisions.KNewOrbit import *

au= 149597871e3
G = 6.67408e-11

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


def ICLVelocitiesGraph(a1, e1, a2, e2, Mstar):
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

    plt.figure()
    plt.subplot(211)
    plt.plot(L/(pi),abs(Rdot1[0,:]), label='rd1 at a')
    plt.plot(L/(pi),RThetadot1[0,:],label='rtd1 at a')
    plt.plot(L/(pi), abs(Rdot2[0, :]), label='rd2 at a')
    plt.plot(L / ( pi), RThetadot2[0, :], label='rtd2 at a')
    plt.legend()
    plt.subplot(212)
    plt.plot(L / ( pi), abs(Rdot1[1, :]), label='rd1 at b')
    plt.plot(L / ( pi), RThetadot1[1, :], label='rtd1 at b')
    plt.plot(L / ( pi), abs(Rdot2[1, :]), label='rd2 at b')
    plt.plot(L / (pi), RThetadot2[1, :], label='rtd2 at b')
    plt.legend()

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
    plt.show()


    return

def ICLVelocities(a1, e1, a2, e2, Mstar):
    N=1000 # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)   #NOTE CHECK L def!! + or -
    R= np.zeros((2,N))  #radial collision points
    C = np.zeros((2, N))   #theta collision points

    for i in range(1,N):
        [[C[0,i], C[1,i]], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    Rdot1=np.zeros(N) # another (2,N) array. Note taken out masses, no differnece!
    Thetadot1=np.zeros(N)
    Rdot2=np.zeros(N)
    Thetadot2=np.zeros(N)

    vrelpar=np.zeros((N))
    for i in range(1,N-1):
        Rdot1[i] = rdot(a1, e1, C[0, i], Mstar)  # another (2,N) array. Note taken out masses, no differnece!
        Thetadot1[i] = thetadot(a1, e1, C[0, i], Mstar)
        Rdot2[i] = rdot(a2, e2, C[0, i] - L[i], Mstar)
        Thetadot2[i] = thetadot(a2, e2, C[0, i] - L[i], Mstar)
        vrelpar[i]=np.sqrt(((Rdot1[i]-Rdot2[i])**2.)+((R[0,i]*(Thetadot1[i]-Thetadot2[i]))**2.))
    #print(L[0]/pi,L[N-1]/pi)
    #print('vrelpar at L=0',vrelpar[0])
    #print('vrelpar at L=2pi', vrelpar[N-1])
    plt.figure()
    plt.semilogy(L[1:N-2]/pi,vrelpar[1:N-2],label='vrelpar')
    plt.xlabel('L/pi')
    plt.ylabel('vrelpar')
    plt.title([a1/au, e1, a2/au, e2])
    plt.show()

    return


ICLVelocities(2 * au, 0.99, 2.1 * au, 0.993, 1.2e30)

#ICLVelocitiesGraph(2 * au, 0.99, 2.1 * au, 0.993, 1.2e30)