from Funcs import *
import matplotlib.pyplot as plt
import numpy as np
from KCollisions.NewOrbit import *
from KCollisions.MinimalOrbitCascade import *
G = 2.982e-27

def erdot(e,R,Thetadot,Mstar): #gives an rdot from r, positive only
    u=G*Mstar
    Rdot=(u/((R**2.)*abs(Thetadot)))*np.sqrt((e**2.)-(((((R**3.)*(Thetadot**2.))/u)-1)**2.))
    return Rdot

def erdotdefdata(e,R,Mstar,n):
    ThetaBottom=(np.sqrt(1-e))*thetakepler(R,Mstar)
    ThetaTop = (np.sqrt(1+e)) * thetakepler(R, Mstar)
    ThetaDot=np.linspace(ThetaBottom,ThetaTop,n+2)
    RDot=erdot(e,R,ThetaDot[1:n+1],Mstar) #ThetaDot[1:n+1]
    eMomentas=np.zeros((2,2*n))
    eMomentas[0,0:n]=ThetaDot[1:n+1]
    eMomentas[1,0:n]=RDot[:]
    eMomentas[:,n:((2*n)+1)]=np.fliplr(eMomentas[:,0:n])
    eMomentas[1,n:((2*n)+1)]=-eMomentas[1,n:((2*n)+1)]
    return eMomentas
def elines(R,Mstar):
    n=10000
    plt.plot(0, R*thetakepler(R, Mstar), 'ro', label='e=0')
    for i in np.arange(0.1, 1.1, 0.1):
        eMomentas = erdotdefdata(i, R, Mstar, n)
        plt.plot(eMomentas[1, :], R * eMomentas[0, :], label='e=%s' % i)
    return

def eGraph(R,Mstar):
    n=10000
    plt.figure()
    elines(R, Mstar)
    plt.xlabel('Rdot')
    plt.ylabel('R*ThetaDot')
    plt.title('Lines of Constant e with R=%s'%R)
    plt.legend()
    plt.show()
    return



def ComGraph(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N):
    K=1

    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    CascadeData=MinimalOrbitCascade(x,a1, e1, s1, a2, e2, s2, Mstar, m1, m2, K, N)

    if CascadeData[N,0,0] > 2:
        LargestGenSize = int(CascadeData[N,0,0])
    else:
        LargestGenSize = 2
    Momenta=np.zeros((N+1,LargestGenSize+1,2)) #keep same format as before, but rdot and thetadot in 3rd dim
    plt.figure()
    for i in range(0,N+1):
        #print ('gen size',CascadeData[i,0,0])
        for j in range (1,int(CascadeData[i,0,0])+1):
            Momenta[i,j,0]=rdot(CascadeData[i,j,3],CascadeData[i,j,4],C-CascadeData[i,j,5],Mstar)  #Check Ls and S!
            Momenta[i, j, 1] = thetadot(CascadeData[i, j, 3], CascadeData[i, j, 4], C - CascadeData[i, j, 5], Mstar)
            plt.plot(Momenta[i,j,0],R*Momenta[i,j,1],'o')#,label='Orbit %s'%CascadeData[i,j,0])
    #plotting efixed
    elines(R,Mstar)

    plt.legend()
    plt.xlabel('R Dot')
    plt.ylabel('R* Theta Dot')
    plt.title('Collision Point %s'%x)
    plt.show()
    return
def SimpleComGraph(a1,e1,s1,a2,e2,s2,Mstar,m1,m2):
    C = np.zeros((2, 1))
    R = np.zeros((2, 1))
    [[C[0], C[1]], [R[0], R[1]]] = CollisionPoints(a1, e1, s1, a2, e2, s2)

    Rdot1=rdot(a1,e1,C-s1,Mstar)
    Thetadot1=thetadot(a1,e1,C-s1,Mstar)
    Rdot2 = rdot(a2, e2, C - s2, Mstar)
    Thetadot2 = thetadot(a2, e2, C - s2, Mstar)

    Rdot3=(m1*Rdot1+m2*Rdot2)/(m1+m2)
    Thetadot3 = (m1 * Thetadot1 + m2 * Thetadot2) / (m1 + m2)
    ThetaKepler=thetakepler(R,Mstar)


    point=['a','b']
    for i in range(0,2):
        #plt.subplot(12(i+1))
        plt.figure()
        plt.plot([0], ThetaKepler[i], 'o', label='e=0')
        plt.plot(Rdot1[i],Thetadot1[i],'o',label='Orbit 1')
        plt.plot(Rdot2[i], Thetadot2[i], 'o', label='Orbit 2')
        plt.plot(Rdot3[i], Thetadot3[i], 'o', label='Orbit 3 ')

        plt.title('CollisionPoint %s' %point[i])
        plt.legend()
        plt.xlabel('Rdot')
        plt.ylabel('Theta Dot')
    plt.show()
    return




ComGraph('a',2,0.99,0,2.1,0.993,2,1.2e30,2e10,2e10,10)
#SimpleComGraph(2*au,0.99,0,2.1*au,0.993,1,1.2e30,2e10,2e10)

#erdotdefdata(0.2,1,1.2e30, 10)
#eGraph(1,1.2e30)