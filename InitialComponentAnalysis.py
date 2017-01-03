import numpy as np
from math import *
from Funcs import *
from NewOrbit import *
import matplotlib.pyplot as plt

au= 1#149597871e3
G = 2.982e-27#6.67408e-11
def InitialComponents(a1,e1,a2,e2,Mstar,m1,m2,N):
    L = np.linspace(0, 2 * pi, N)   #testing
    R= np.zeros((2,N))  #radial collision points
    C = np.zeros((2, N))   #theta collision points

    for i in range(1,N):
        [[C[0,i], C[1,i]], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, -L[i])


    #a data stuff
    ThetaKepler=np.sqrt((G*Mstar)/(np.power((R[:,:]),3)))
    Rdot1=rdot(a1,e1,C[:,:],m1+Mstar)
    ThetadotK1=thetadot(a1,e1,C[:,:],m1+Mstar)-ThetaKepler
    Rdot2 = rdot(a2, e2, C[:, :]-L, m2 + Mstar)
    ThetadotK2 = thetadot(a2, e2, C[:, :]-L, m2 + Mstar)-ThetaKepler
    abtitle=['a','b']
    for i in range(0,2): #loop plotting for both collision points
    #rdot
        plt.figure(1+i)
        plt.subplot(221)
        plt.plot(L[1:N]/pi, Rdot1[i,1:N], label="Rdot1")
        plt.plot(L[1:N]/pi, Rdot2[i,1:N], label="Rdot2")
        plt.legend() #attaching a legend to the graph
        plt.title('Rdot against Lambda/pi')
        plt.xlabel('Lambda/pi, Sepetation of the Orbits')
        plt.ylabel('Rdot')
        plt.xlim([0,2])
        plt.subplot(222)
        plt.plot(L[1:N]/pi, np.multiply(Rdot1[i,1:N],Rdot2[i,1:N]), label="Rdot1*Rdot2")
        plt.xlim([0,2])
        #athetadot
        plt.subplot(223)
        plt.plot(L[1:N]/pi, ThetadotK1[i,1:N], label="Thetadotk1")
        plt.plot(L[1:N]/pi, ThetadotK2[i,1:N], label="Thetadotl2")
        plt.legend()  # attaching a legend to the graph
        plt.title('Thetadotk against Lambda/pi')
        plt.xlabel('Thetadot')
        plt.xlim([0, 2 ])
        plt.subplot(224)
        plt.plot(L[1:N]/pi, np.multiply(ThetadotK1[i,1:N], ThetadotK2[i,1:N]), label="Rdot1*Rdot2")
        plt.xlim([0,2])
        plt.suptitle('collision point %s'%(abtitle[i]))

    #graph of the collision points

    plt.figure(3)
    plt.plot(L[1:N] / pi, R[0, 1:N], label="Collision RAdius a")
    plt.plot(L[1:N] / pi, R[1, 1:N], label="Collision Radius b")
    plt.plot([0,2],[1,1], label="1au")
    plt.legend()
    plt.title('Collision Points')
    plt.show()

    return

InitialComponents(2*au,0.99,2.1*au,0.993,1.2e30,2e10,2e10,100)




