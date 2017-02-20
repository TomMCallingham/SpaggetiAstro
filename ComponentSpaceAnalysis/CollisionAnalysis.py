import numpy as np
from numpy import linalg as nplg
import matplotlib.pyplot as plt
from KCollisions.NewOrbit import *
Mstar=1.2e30
au=1.496e11
year=3.154e+7
def CollisionAngle(x,a1,e1,s1,a2,e2,s2):
    N=1000
    L = s1 - s2
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    rd1=rdot(a1,e1,C-s1,Mstar)
    td1 = thetadot(a1, e1, C - s1, Mstar)
    V1=np.array([rd1,R*td1])
    v1=nplg.norm(V1)

    rd2 = rdot(a2, e2, C - s2, Mstar)
    td2 = thetadot(a2, e2, C - s2, Mstar)
    V2 = np.array([rd2, R * td2])
    v2=nplg.norm(V2)

    Vrel=V2-V1
    Momentas=np.zeros((3,N))
    H=np.linspace(0,1,N)


    Angles1=np.zeros(N)
    Angles2=np.zeros(N)

    for i in range(0,N):
        Momentas[0:2, i] = V1 + H[i] * Vrel
        Momentas[2,i]=nplg.norm(Momentas[0:2,i])
        Angles1[i]=np.arccos(np.dot(Momentas[0:2,i],V1)/(v1*Momentas[2,i]))
        Angles2[i] = np.arccos(np.dot(Momentas[0:2, i], V2) / (v2 * Momentas[2,i]))
    '''
    plt.figure()
    plt.plot(H,)
    '''

    plt.figure()
    plt.plot(H,Angles1,label='angles1')
    plt.plot(H,Angles2, label='Angles2')
    '''
    #trial
    plt.plot([0,1],[0,Angles1[N-1]], label='guess')
    plt.plot([0, 1], [Angles2[N - 1], 0], label='guess')
    plt.legend()

    plt.figure()
    plt.plot(H[1:N],1/np.sin(Angles1[1:N]), label='cross volume?')
    plt.plot(H[1:N], 1 / np.sin(Angles2[1:N]), label='cross volume?')
    plt.legend()
    plt.title('L=%s'%L)
    '''
    plt.show()

    return

def oldCollisionAngle(x,a1,e1,s1,a2,e2,s2):
    N=1000
    L = s1 - s2
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    rd1=rdot(a1,e1,C-s1,Mstar)
    td1 = thetadot(a1, e1, C - s1, Mstar)
    V1=np.array([[rd1],[R*td1]])
    v1=nplg.norm(V1)

    rd2 = rdot(a2, e2, C - s2, Mstar)
    td2 = thetadot(a2, e2, C - s2, Mstar)
    V2 = np.array([[rd2], [R * td2]])
    v2=nplg.norm(V2)

    Vrel=V2-V1
    Momentas=np.zeros((3,N))
    H=np.linspace(0,1,N)

    Momentas[0:2,:]=V1 +H*Vrel
    Angles1=np.zeros(N)
    Angles2=np.zeros(N)

    for i in range(0,N):
        Momentas[2,i]=nplg.norm(Momentas[0:2,i])
        Angles1[i]=np.arccos(np.dot(Momentas[0:2,i],V1)/(v1*Momentas[2,i]))
        Angles2[i] = np.arccos(np.dot(Momentas[0:2, i], V2) / (v2 * Momentas[2,i]))
    '''
    plt.figure()
    plt.plot(H,)
    '''

    plt.figure()
    plt.plot(H,Angles1,label='angles1')
    plt.plot(H,Angles2, label='Angles2')
    '''
    #trial
    plt.plot([0,1],[0,Angles1[N-1]], label='guess')
    plt.plot([0, 1], [Angles2[N - 1], 0], label='guess')
    plt.legend()

    plt.figure()
    plt.plot(H[1:N],1/np.sin(Angles1[1:N]), label='cross volume?')
    plt.plot(H[1:N], 1 / np.sin(Angles2[1:N]), label='cross volume?')
    plt.legend()
    plt.title('L=%s'%L)
    '''
    plt.show()

    return



def MomentaData(x,a1,e1,s1,a2,e2,s2):
    dfSize=1#1/(1.496e11) #1 meter
    Rbeam=1000*dfSize

    N = 1000
    L = s1 - s2
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]


    #setting up the parent data
    rd1 = rdot(a1, e1, C - s1, Mstar)
    td1 = thetadot(a1, e1, C - s1, Mstar)
    V1 = np.array([rd1, R * td1])
    v1 = nplg.norm(V1)
    rd2 = rdot(a2, e2, C - s2, Mstar)
    td2 = thetadot(a2, e2, C - s2, Mstar)
    V2 = np.array([rd2, R * td2])
    v2 = nplg.norm(V2)
    Vrel = V2 - V1
    vrel=nplg.norm(Vrel)

    #Creating the Momentas : rd,R*td,v,vrel1,vrel2
    MomentasData = np.zeros((11, N))#,dtype=float) # 0-rd,1-R*td,2-v,3-vrel1,4-vrel2,5-angle1,6-angle2,7-e,8-TimePeriod,9-T1,10-T2
    H = np.linspace(0, 1, N)



    #Finding The Collision Angles, skipping parents
    for i in range(0,N):
        MomentasData[0:2, i] = V1 + H[i] * Vrel  # rd,td
        MomentasData[2, i] = nplg.norm(MomentasData[0:2, i])  # v
        MomentasData[3, i] = H[i]*vrel #nplg.norm(MomentasData[0:2, i] - V1)  # vrel1
        MomentasData[4, i] = (1-H[i])*vrel#nplg.norm(MomentasData[0:2, i] - V2)  # vrel2
        MomentasData[5,i]=np.arccos(np.vdot(MomentasData[0:2,i],V1)/(v1*MomentasData[2,i])) #angle 1
        MomentasData[6, i] = np.arccos(np.vdot(MomentasData[0:2, i], V2) / (v2 * MomentasData[2,i])) #angle2
    # e calc
    MomentasData[7,:]=newe(MomentasData[0,:],MomentasData[1,:]/R,R,Mstar) #eccentricity
    #calculating the Time period
    MomentasData[8,:]=((2*pi*(R**6.))/((G*Mstar)**2.))*((((MomentasData[1,:]/R)**3.))/np.sqrt(((1-(MomentasData[7,:]**2.))**(3.))))
    #calculating collision time periods
    MomentasData[9,1:N-1]=(pi*Rbeam)/(2*MomentasData[2,1:N-1])
    MomentasData[10,1:(N-1)]=MomentasData[9,1:N-1]/ abs(np.sin(MomentasData[6, 1:N-1])) #T2
    MomentasData[9, 1:N-1] = MomentasData[9, 1:N-1] / abs(np.sin(MomentasData[5, 1:N-1]))#T1

    alpha1=abs(np.arccos(np.vdot(V1,Vrel)/(v1*vrel)))
    alpha2 = abs(np.arccos(np.vdot(V2, Vrel) / (v2 * vrel)))
    '''
    #Eccentricity Graph
    plt.figure()
    plt.plot(H,MomentasData[7,:], label='eccentricity')
    plt.legend()

    #TimePeriodGraph
    plt.figure()
    plt.plot(H, MomentasData[8, :]/year, label='TimePeriod in years')
    #plt.plot(np.array([0,1]),((2*pi)/sqrt(G*Mstar))*np.array([a1**1.5,a2**1.5]),'o',label='check')
    plt.legend()

    #Collision Time Graph
    plt.figure()
    plt.plot(H[1:N-1], MomentasData[9, 1:(N-1)]/year, label='T1 in year')
    plt.plot(H[1:N-1],MomentasData[10,1:(N-1)]/year, label='T2 in years')
    plt.legend()

    #Fractional Time Graph
    plt.figure()
    plt.plot(H[1:N - 1], MomentasData[9, 1:(N - 1)]/MomentasData[8,1:(N-1)], label='T1/T')
    plt.plot(H[1:N - 1], MomentasData[10, 1:(N - 1)]/MomentasData[8,1:(N-1)], label='T2/T')
    plt.legend()

    #vplot
    plt.figure()
    plt.plot(H,MomentasData[2,:],label='v?')
    plt.legend()

    #Anlgeplot
    plt.figure()
    plt.plot(H[1:N-1],MomentasData[5, 1:N - 1]/pi,label='Angle1/pi')
    plt.plot(H[1:N - 1], MomentasData[6, 1:N - 1]/pi, label='Angle2/pi')
    plt.legend()
    #singraph
    plt.figure()
    plt.plot(H[1:N - 1], np.sin(MomentasData[5, 1:N - 1]), label='sin angle1')
    plt.plot(H[1:N - 1], np.sin(MomentasData[6, 1:N - 1]), label='sin angle2')
    plt.legend()

    #1/singraph
    plt.figure()
    plt.plot(H[1:N-1],1/abs(np.sin(MomentasData[5, 1:N - 1])), label='1/sin angle1')
    plt.plot(H[1:N - 1], 1/abs(np.sin(MomentasData[6, 1:N - 1])), label='1/sin angle2')
    plt.legend()

    # vrel
    plt.figure()
    plt.plot(H[1:N - 1], MomentasData[3, 1:N - 1] , label='vrel1')
    plt.plot(H[1:N - 1], MomentasData[4, 1:N - 1] , label='vrel2')
    plt.legend()

    #vrel/sin
    plt.figure()
    plt.plot(H[1:N - 1], MomentasData[3, 1:N - 1] / abs(np.sin(MomentasData[5, 1:N - 1])), label='vrel1/sin angle1')
    plt.plot(H[1:N - 1], MomentasData[4, 1:N - 1] / abs(np.sin(MomentasData[6, 1:N - 1])), label='vrel 2/sin angle2')
    plt.legend()
    '''
    #vrel/vsin
    plt.figure()
    plt.plot(H[1:N - 1], MomentasData[3, 1:N - 1] / (MomentasData[2, 1:N - 1] *abs(np.sin(MomentasData[5, 1:N - 1]))), label='vrel1/vsin angle1')
    plt.plot(H[1:N - 1], MomentasData[4, 1:N - 1] / (MomentasData[2, 1:N - 1] *abs(np.sin(MomentasData[6, 1:N - 1]))), label='vrel 2/vsin angle2')
    plt.plot(np.array([0.5,1]),(1/sin(alpha1))*np.array([1,1]),label='sincheck')
    plt.legend()


    #Rcol
    Rcol=np.zeros((2,N))
    MomentaCol=np.zeros((2,N)) #orbit dependand part, 1 for each parent
    PConstant=(3 * ((G * Mstar) ** 4.) * pi)/(256*Rbeam*(R**12.))
    ParentConstant=np.zeros((2))
    ParentConstant[0]=1#/(v1*((td1**3.))) #parent1 constant
    ParentConstant[1] = 1# / (v2 * ((td2 ** 3.)))  # parent1 constant
    rsum=20 #largest fragments

    MomentaCol[0,:]=(MomentasData[3,:]/(MomentasData[2,:]*((MomentasData[1,:]/R)**3.)*abs(np.sin(MomentasData[5,:])))) #parent1
    MomentaCol[1, :] = (MomentasData[4, :] / (MomentasData[2, :] * ((MomentasData[1, :] / R) ** 3.) * abs(np.sin(MomentasData[6, :]))))#paernt2

    #year=
    Rcol[0,:]=PConstant*ParentConstant[0]*MomentaCol[0,:]*(dfSize ** 2.) * (rsum**2.)
    Rcol[1, :] = PConstant * ParentConstant[1] * MomentaCol[1, :] * (dfSize ** 2.) *( rsum**2.)

    plt.figure()
    plt.plot(H,np.log10(Rcol[0,:]*year),label='log10 (Rcol year) Parent1')
    plt.plot(H, np.log10(Rcol[1, :]*year), label='log10(Rcoll yaer) Parent2')
    plt.legend()

    plt.show()

    return MomentasData





    #PConstant = (3 * ((G * Mstar) ** 4.) * pi) / 256 * Rbeam * (R ** 12.)
    #ParentProb=np.zeros(2) # 0 index is parent 1, 1 index parent 2
    #ParentProb[0]=1/(v1*(td1**3.))
    #ParentProb[1] = 1 / (v2 * (td2 ** 3.))"""



MomentaData('b',2*au,0.99,0,2.1*au,0.993,2.1)
#CollisionAngle('a',2,0.99,0,2.1,0.993,1.5)