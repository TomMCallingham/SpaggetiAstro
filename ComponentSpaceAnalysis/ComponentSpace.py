from Funcs import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from KCollisions.NewOrbit import *
from KCollisions.MinimalOrbitCascade import *
from ComponentSpaceAnalysis.EccentricitySolver import *
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
    K=0.9
    Mstar=1.2e30
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
    colours = cm.rainbow(np.linspace(0, 1, N+1))

    for i in range(0,N+1): #generation loop
        #print ('gen size',CascadeData[i,0,0])
        for j in range (1,int(CascadeData[i,0,0])+1):
            Momenta[i,j,0]=rdot(CascadeData[i,j,3],CascadeData[i,j,4],C-CascadeData[i,j,5],Mstar)  #Check Ls and S!
            Momenta[i, j, 1] = thetadot(CascadeData[i, j, 3], CascadeData[i, j, 4], C - CascadeData[i, j, 5], Mstar)
            #plt.plot(Momenta[i,j,0],R*Momenta[i,j,1],'o' )#,color=colours[i])#,label='Orbit %s'%CascadeData[i,j,0])
    #plotting efixed
        plt.plot(Momenta[i,1:int(CascadeData[i,0,0])+1,0],R*Momenta[i,1:int(CascadeData[i,0,0])+1,1],'o',color=colours[i], label='gen %s'%i )#,color=colours[i])#,label='Orbit %s'%CascadeData[i,j,0])
    elines(R,Mstar)

    [rdmin, tdmin]=eSolver(Momenta[0,1,0],Momenta[0,1,1],Momenta[0,2,0],Momenta[0,2,1],thetakepler(R,Mstar),R)
    plt.plot(rdmin,R*tdmin,'ks', label='Best')
    #emin=newe(rdmin,tdmin,R,Mstar)
    #print('emincheck',emin)

    plt.legend()
    plt.xlabel('R Dot')
    plt.ylabel('R* Theta Dot')
    plt.title('Collision Point %s'%x)
    plt.show()
    return

def MineComGraph(x,a1,e1,s1,a2,e2,s2,Mstar):
    K=1
    Mstar=1.2e30
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    rd1=rdot(a1,e1,C-s1,Mstar)
    td1=thetadot(a1,e1,C-s1,Mstar)
    rd2 = rdot(a2, e2, C - s2, Mstar)
    td2 = thetadot(a2, e2, C - s2, Mstar)
    tdk=thetakepler(R, Mstar)

    plt.figure()
    [rdmin, tdmin] = eSolver(rd1, td1, rd2, td2, tdk, R)
    emin=newe(rdmin,tdmin,R, Mstar)
    plt.plot(rdmin, R * tdmin, 'rs', label='Best')
    elines(R, Mstar)
    plt.plot(rd1, R*td1, 'o', label='Orbit 1')
    plt.plot(rd2, R*td2, 'o', label='Orbit 2')
    plt.plot(np.array([rd1,rd2]),R*np.array([td1,td2]),'b')


    plt.title('CollisionPoint %s, emin=%s' %(x,emin))
    plt.legend()
    plt.xlabel('Rdot')
    plt.ylabel('Theta Dot')
    plt.show()
    return





ComGraph('a',2,0.99,0,2.1,0.993,1.8,1.2e30,2e10,2e10,10)
#MineComGraph('a',2*au,0.99,0,2.1*au,0.993,1,1.2e30)

#erdotdefdata(0.2,1,1.2e30, 10)
#eGraph(1,1.2e30)

#esolvertest('a',2,0.99,0,2.1,0.993,2)
