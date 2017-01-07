import numpy as np
from Funcs import *
from KCollisions.NewOrbit import *
import matplotlib.pyplot as plt

def KEnergy(a1,e1,s1,a2,e2,s2,Mstar,m1,m2):
    L=s1-s2
    N=10
    K=np.linspace(0.1,1,N)
    [C, R] = CollisionPoints(a1, e1, 0, a2, e2, -L)
    Rdot1=np.zeros((2,1))
    Rdot2 = np.zeros((2, 1))
    Thetadot1 = np.zeros((2, 1))
    Thetadot2 = np.zeros((2, 1))
    Rdot3=np.zeros((2,N))
    Thetadot3 = np.zeros((2, N))
    E = np.zeros((2, N))
    for j in range(0,2): #2 Collision Points
        for i in range(0,N):
            [a3, e3, s3, m3]=NewOrbit(a1,e1,0,a2,e2,-L,Mstar,m1,m2,K[i],C[j],R[j])
            Rdot3[j,i]=rdot(a3,e3,s3-L,m3) #at colision as well? AT here
            Thetadot3[j,i]=thetadot(a3,e3,s3-L,m3)
            #E[j,i]=m3*SEnergy(Rdot3[j, i], Thetadot3[j, i], R[j], Mstar)-(((1-K[i])*(m1+m2)*G*Mstar)/R[j])
            E[j,i]=(m3*((Rdot3[j,i]**2.)+((R[j]*Thetadot3[j,i])**2.)))/2

        Rdot1[j] = rdot(a1, e1, s1 - L, m1)
        Thetadot1[j] = thetadot(a1, e1, s1 - L, m1)
        Rdot2[j]=rdot(a2, e2, s2 - L, m2)
        Thetadot2[j] = thetadot(a2, e2, s2 - L, m2)



    # Calculating Ebefore
    #Ebefore= m1 * SEnergy(Rdot1, Thetadot1, R, Mstar) + m2 * SEnergy(Rdot2, Thetadot2, R, Mstar)
    Ebefore=((m1*((Rdot1**2.)+((R*Thetadot1)**2.)))+(m2*((Rdot2**2.)+((R*Thetadot2)**2.))))/2
    print('Ebefore'), print(Ebefore)
    print('Rdot3'), print(Rdot3)


    plt.figure()
    plt.subplot(121)
    plt.plot(K,E[0,:],label='Energy After at a')
    plt.plot([0,1],Ebefore[:,0],label='Energy Before at a')
    plt.legend()

    plt.subplot(122)
    plt.plot(K, E[1, :], label='Energy After at b')
    plt.plot([0, 1], Ebefore[:,1], label='Energy Before at b')
    plt.legend()

    plt.show()
    return

KEnergy(2,0.99,0,2.1,0.993,2,1.2e30,2e10,2e10)

