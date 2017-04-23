import numpy as np
import matplotlib.pyplot as plt
from numba import jit
from Funcs import *
from ComponentSpaceAnalysis.EccentricitySolver import *
from KCollisions.KNewOrbit import *
#Units

au = 1.496e11
year = 3.154e+7
Mstar = 1.2e30
Lsun=3.828e26
Msun=1.98e30
WdRad=8.75e6
RocheRad=0.013*au#2.5e9
Q=1
#m=1, A=1 Assume Spherical
Density=2e3#1.82e3 #from paper
c=3e8
#Rfrag=10**(-1.)
'''
qmin=5.86e-5
qmax=0.013
qini=((qmax+(19*qmin))/20)*au
aini=10*au
eini=1-(qini/aini)
print('eini',eini)
'''
aini=10*au
eini=0.9993
Tini=5e6

econstant = ((5 * Q) / (8 * np.pi * (c ** 2.))) * (3 / Density)
aconstant = (Q / (4 * np.pi * (c ** 2.))) * (3 / Density)

def Luminosity(t): #note time in years
    L=3.26*Lsun*((0.1+t*(10**-6.))**(-1.18))
    return L
def GraphLuminosity():
    Time=np.linspace(0,9*(10**9.),1e7)
    L=Luminosity(Time)
    plt.figure()
    plt.loglog(Time,L/Lsun,label='Luminosity in Lsun')
    plt.ylabel('Luminosity (Lsun)')
    plt.xlabel('Time')
    plt.title('WD Luminosity against time')
    plt.legend()
    plt.show()
    return

@jit
def PREvolve(a,e,Tyears,Tini,Rfrag):
    Timestep=10**1.  #years
    #Tyears=10**8.
    T=int(Tyears/Timestep)
    #print('T=',T)
    Timestep=Timestep*year

    Data=np.zeros((3,T+1)) #a,e one left blank!
    Data[0,:]=np.arange(0, T+1)*Timestep
    Data[1:,0]=[a,e]#,a*(1-e)]
    L=Luminosity((np.arange(0,T+1))+Tini)
    #plt.plot(Data[0, :], L, label='L')

    Tmini=T*Timestep
    result='unfinished'
    for t in range(0, T):
        Data[1, t + 1] = Data[1, t] - Timestep * (
            (aconstant/Rfrag) * L[t] * (2 + (3 * (Data[2, t] ** 2.)))) / (
                                                  Data[1, t] * ((1 - (Data[2, t] ** 2.)) ** 1.5))  # a
        Data[2, t + 1] = Data[2, t] - Timestep * (((econstant/Rfrag) * L[t ] * Data[2, t]) / (
        (Data[1, t] ** 2.) * (np.sqrt(1 - (Data[2, t] ** 2.)))))  # e
        #Data[3,t+1]=Data[1,t+1]*(1-Data[2,t+1])

        if Data[2, t + 1] < 0: #shrunk e<0.2
            Tmini=t*Timestep
            result='shrunk'
            Data[2, t + 1]=0
            Data[1,t+2:]=Data[1, t + 1]
            Data[2, t+2:] = Data[2, t + 1]
            break
        '''
        if Data[1, t + 1]*(1+Data[2,t+1]) < RocheRad/100: #awithin the Roche Rad Should be a*(1+e), but nearly the same thing
            #Data[1, t + 1] = RocheRad
            #Data[2, t + 1] = 0
            #Data[3,t+1]=
            Tmini=t*Timestep
            result='within Roche RAd'
            Data[1, t+2:] = Data[1, t + 1]
            Data[2, t+2:] = Data[2, t + 1]
            #Data[3,t+2:]=Data[3,t+1]
            break
        '''



    print(result)#,Tmini)

    return (Data,Tmini)


def GraphPrEvolve(a,e,T,Tini):
    plt.figure(1) #a
    plt.figure(2) #e
    for i in range(-3,3):
        Rfrag=10**i
        (Data, Tmini)=PREvolve(a,e,T,Tini,Rfrag)
        plt.figure(1)
        plt.loglog(Data[0,:]/year,Data[1,:]/au,label='Rfrag=10^%s'%i) #a plot
        plt.figure(2)
        plt.semilogx(Data[0, :]/year, Data[2, :], label='Rfrag=10^%s'%i)  #e plot
    plt.legend()
    #plt.title(' Radiation shrinking Eccentricity, Tini=%.E'%Tini)
    plt.ylabel('eccentricity')
    plt.xlabel('Time/years')
    plt.xlim([1e4, 1e8])
    plt.title('a=%s au, e=%s'%(a/au,e))
    plt.figure(1)
    #plt.title(' Radiation shrinking a Tini=%.E'%Tini)
    plt.ylabel('a/au')
    plt.xlabel('Time/years')
    plt.title('a=%s au, e=%s' % (a/au, e))
    plt.ylim([(RocheRad/au)*(1/10),(2*a)/au])
    plt.xlim([1e4,1e8])
    plt.legend()
    plt.show()
    return

def GraphPREvolveVar(R,rd,td,T,Tini):
    e=newe(rd,td,R,Mstar)
    a=newa(rd,td,R,Mstar)
    GraphPrEvolve(a, e, T, Tini)
    return

def eminPR(a1,e1,a2,e2,L,x):
    CollisionData = CollisionPoints(a1, e1, 0, a2, e2, L)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    rd1 = rdot(a1, e1, C , Mstar)
    td1 = thetadot(a1, e1, C , Mstar)
    rd2 = rdot(a2, e2, C - L, Mstar)
    td2 = thetadot(a2, e2, C - L, Mstar)
    [rdmin, tdmin] = eSolver(rd1, td1, rd2, td2, R)
    T=10**9.
    Tini=0.5e9
    GraphPREvolveVar(R,rdmin,tdmin,T,Tini)
    return

@jit
def GraphTimePREvolve(a1, e1, a2, e2,x,Rfrag):
    j=1
    if x=='a':
        j=0


    MineData = LambdaeData(a1, e1, a2, e2)
    N = np.size(MineData[0, 0, :])
    L = np.linspace(0, 2 * pi, N + 2)  # NOTE CHECK L def!! + or -
    L = L[1:N + 1]
    R = np.zeros((2, N))  # radial collision points
    C = np.zeros((2, N))  # theta collision points

    for i in range(0, N):
        [[C[0, i], C[1, i]], [R[0, i], R[1, i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])
    Tmini=np.zeros((N))
    amat=np.zeros((N))
    amat=newa(MineData[j,1,:],MineData[j,2,:],R[j,:],Mstar)
    for i in range(0,N):
        #print('e',MineData[j,2,i])
        #print('a',amat[i])

        Tmini[i]=PREvolve(amat[i], MineData[j,0,i], 10**8., Tini, Rfrag)[1]

    print(Tmini[0],Tmini[N-1])
    plt.figure()
    plt.plot(L/pi,Tmini/year)
    plt.xlabel('Lambda/pi')
    plt.ylabel('Time taken to e<0.2 or hit WD ')
    plt.show()
    return
#GraphPrEvolve(aini,eini,10**8,0.5e9)
#GraphTimePREvolve(2*au,0.99,2.1*au,0.993,'a',10)
#GraphLuminosity()
#eminPR(2*au,0.99,2.1*au,0.993,1.5,'a')


#print(0.005*au)

def changeperorbit(tini, a1, e1):#,a2,e2):

    Tp=(2*pi)*np.sqrt((a1 ** 3.) / (G * Mstar))

    Rda =- Tp * (
        (aconstant ) * Luminosity(tini) * (2 + (3 * (e1 ** 2.)))) / (
             a1 * ((1 - (e1 ** 2.)) ** 1.5))  # R*da
    Rde =- Tp * ((econstant * Luminosity(tini) * e1) / (
        (a1 ** 2.) * (np.sqrt(1 - (e1 ** 2.))))) # R*de

    print('Rda=',Rda)
    print('Rde=',Rde)

    return



def ringdispresion(a,e,f,Tyears,Tini):

    #plt.figure(1) #change in a
    #plt.figure(2) # change in e
    plt.figure(3) #pericentre
    for i in range(-3,3):
        Rfrag=10**i
        (Data, Tmini)=PREvolve(a,e,Tyears,Tini,Rfrag)
        T=np.size(Data[0,:])-1
        #print(T)
        DataAnalysis=np.zeros((1,T+1))
        for t in range(0,T+1):
            DataAnalysis[0,t]=npr(Data[1,t],Data[2,t],f) #(1-Data[2,t])*Data[1,t]#
        '''
        plt.figure(1)
        plt.loglog(Data[0, :] / year, (a-Data[1,:]) , label='Rfrag=10^%s' % i)  # change in a plot
        plt.figure(2)
        plt.loglog(Data[0, :] / year, (e-Data[2,:]), label='Rfrag=10^%s' % i)  # change in a plot
        '''
        plt.figure(3)
        plt.loglog(Data[0, :] / year, (DataAnalysis[0, 0]-DataAnalysis[0,:]), label='Rfrag=10^%s' % i)  # change in q


    # plt.title(' Radiation shrinking Eccentricity, Tini=%.E'%Tini)
    plt.ylabel('change in q /m')
    plt.xlabel('Time/years')
    plt.title('a=%s au, e=%s, Tini=%.f' % (a / au, e,Tini))
    plt.legend()
    '''
    plt.figure(1)
    # plt.title(' Radiation shrinking a Tini=%.E'%Tini)
    plt.ylabel('change in a / m')
    plt.xlabel('Time/years')
    plt.title('a=%s au, e=%s, Tini=%.f' % (a / au, e, Tini))
    plt.legend()
    plt.figure(2)
    # plt.title(' Radiation shrinking a Tini=%.E'%Tini)
    plt.ylabel('change in e')
    plt.xlabel('Time/years')
    plt.title('a=%s au, e=%s, Tini=%.f' % (a / au, e, Tini))
    plt.legend()
    '''

   # plt.loglog([1,Tyears],[10**3.,10**3.], label='Km')


    plt.show()
    print('finished graphing')
    return

#changeperorbit(5e9,2.5*au,0.995)
#ringdispresion(2.1*au,0.995,10**7.,0.5e9)
ringdispresion(2.5 * au, 0.998,0, 10 ** 8., 0.5e9)
#GraphPrEvolve(2.5*au,0.998,10**8,0.5e9)



