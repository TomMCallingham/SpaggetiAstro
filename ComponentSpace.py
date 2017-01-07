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
    ThetaDot=np.linspace(ThetaBottom,ThetaTop,n)
    RDot=erdot(e,R,ThetaDot,Mstar) #ThetaDot[1:n+1]
    eMomentas=np.zeros((2,n))
    eMomentas[0,:]=ThetaDot[:]
    eMomentas[1,:]=RDot[:]
    #print(eMomentas)
    return eMomentas




def ComGraph(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N):
    K=1
    C=np.zeros((2,1))
    R=np.zeros((2,1))
    [[C[0], C[1]], [R[0], R[1]]] = CollisionPoints(a1, e1, s1, a2, e2, s2) #SIGNS MUST BE WRONG THO

    CascadeData=MinimalOrbitCascade(a1, e1, s1, a2, e2, s2, Mstar, m1, m2, K, N)


    plt.figure()
    if CascadeData[N,0,0] > 2:
        LargestGenSize = int(CascadeData[N,0,0])
    else:
        LargestGenSize = 2
    Momenta=np.zeros((N+1,LargestGenSize+1,2)) #keep same format as before, but rdot and thetadot in 3rd dim
    for i in range(0,N+1):
        print ('gen size',CascadeData[i,0,0])
        for j in range (1,int(CascadeData[i,0,0])+1):
            Momenta[i,j,0]=rdot(CascadeData[i,j,3],CascadeData[i,j,4],C[0]-CascadeData[i,j,5],Mstar)  #Check Ls and S!
            Momenta[i, j, 1] = thetadot(CascadeData[i, j, 3], CascadeData[i, j, 4], C[0] - CascadeData[i, j, 5], Mstar)
            plt.plot(Momenta[i,j,0],Momenta[i,j,1],'o',label='Orbit %s'%CascadeData[i,j,0])

    #plotting efixed
    e0=0.2
    n=5
    efixeda=erdotdefdata(e0,R[0],Mstar, n)
   # print(efixeda)
    #efixedb = erdotdefdata(e0, R[1], Mstar, n)
    plt.plot([efixeda[0,:]],[efixeda[1,:]],label='e=0.2')




    #plt.legend()
    plt.xlabel('R Dot')
    plt.ylabel('Theta Dot')


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




ComGraph(2,0.99,0,2.1,0.993,3,1.2e30,2e10,2e10,4)
#SimpleComGraph(2*au,0.99,0,2.1*au,0.993,1,1.2e30,2e10,2e10)

#erdotdefdata(0.2,1,1.2e30, 10)
