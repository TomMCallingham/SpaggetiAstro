from Evolution.RcolTContact import *
from Functions import *
from numba import jit, int32, float64

rmin = 1
rmax = 100
TotalVol = 1e15  # ((4*pi)/3)*(Rast**3.)
# Values
a1 = 2.5 * au
e1 = 0.998
a2 = 2.4 * au
e2 = 0.997
I = 1 * ((2 * pi) / 360)
Rbeam = 1000

torstep=1e-1
Tormax=1.5e4
Torload='TorDist1m.npy'

@jit
def InitialDistCalc(J):
    IniDist = np.zeros((2, J + 1))  # number,rad
    IniDist[1, 1:J + 1] = np.linspace(rmin, rmax, J)  # rad
    for f in range(1, J + 1):
        IniDist[0, f] = (IniDist[1, f]) ** (-3.5)  # Number in size range
        IniDist[0, 0] += IniDist[0, f] * ((IniDist[1, f]) ** 3.) * ((4 * pi) / 3)  # Total Volume
    # print('TotalVol unscaled',Dist[0,0,0])
    IniDist[0, :] = ((TotalVol / IniDist[0, 0])) * IniDist[0, :]
    return IniDist


@jit
def TorStaticEvolution( J):
    J=int(J)
    # Setup
    # Data Matrix Setup
    Tornumber = int(Tormax / torstep)
    TorDist = np.zeros((3, J + 1, Tornumber + 1))  # P1P2C,Frag,Tor
    IniDist = InitialDistCalc(J)
    # print('IniDist',IniDist[0, :])
    TorDist[0, :, 0] = IniDist[0, :]
    TorDist[1, :, 0] = IniDist[0, :]
    # CollisionRate Setup
    CrossSecMat = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            CrossSecMat[i, j] = (IniDist[1, i] + IniDist[1, j]) ** 2.

    # Loop
    print('Evolving...')
    for t in range(0, Tornumber):  # Time Loop
        TorDist[:, 1:J + 1, t + 1] = TorDist[:, 1:J + 1, t]

        for f1 in range(1, J + 1):
            for f2 in range(1, J + 1):  # note mirrored problem!
                colnumber = torstep * CrossSecMat[f1, f2] * TorDist[0, f1, t] * TorDist[1, f2, t]*1e-17  # Calc the number
                # print(colnumber)
                # Remove from parents
                TorDist[0, f1, t + 1] -= colnumber
                TorDist[1, f2, t + 1] -= colnumber
                # add to Child
                TorDist[2, f1, t + 1] += colnumber
                TorDist[2, f2, t + 1] += colnumber

    for f in range(1, J + 1):
        TorDist[:, 0, 1:] += TorDist[:, f, 1:] * ((4 * pi) / 3) * (IniDist[1, f] ** 3.)
    print('Evolved')
    return TorDist

def SaveTorDist(  J):
    TorDist=TorStaticEvolution(  J)
    print('saving')
    np.save('TorDist.npy', TorDist)#,CrossDataminus)
    print('saved')
    return


def TorMassChange():
    loaddist='TorDist.npy'
    print('Loading Data...')
    Dist = np.load(Torload)  # P1P2C,Frag,Tor
    print('Loaded')
    Tnumber = np.size(Dist[0,0,:])
    Tor=np.arange(0,Tormax+torstep,torstep)
    print('Plotting...')
    '''
    plt.figure()
    plt.plot(Tor,Dist[0,0,:], label='P1 Vol')
    plt.plot(Tor, Dist[1, 0, :], label='P2 Vol')
    plt.plot(Tor, Dist[2, 0, :], label='C Vol')
    plt.legend()
    plt.xlabel('Tor')
    plt.ylabel('Vol')
    '''
    #Averaged Evolution Graph

    Lambav=lambdatimeaverage(1e4)*1e17
    plt.figure()
    plt.plot(Tor/Lambav, Dist[2, 0, :]/(2*TotalVol), label='C Vol')
    plt.legend()
    plt.xlabel('Time yrs')
    plt.ylabel('Vol')
    #plt.xlim([0,1e8])


    #Sizes P1 Graph
    '''
    IniDist = InitialDistCalc(J)
    plt.figure()
    for f in range(1,J+1):
        plt.semilogx(Tor,Dist[0,f,:]*(IniDist[1,f]**3.5),label='r=%s'%IniDist[1,f])
    plt.title('P1 Fragments over time')
    plt.legend()
    plt.xlabel('Tor')
    plt.ylabel('Vol')
    '''
    print('finshed')
    return

def PrecessingEvolutionGraphs():
    print('Loading Data...')
    CrossData = np.load('CrossData.npy') #t, lambda=s1-s2,s1,s2,+/-1,a=0 b=1
    TorDist=np.load(Torload)  # P1P2C,Frag,Tor
    print('Loaded')
    M=np.size(CrossData[0,:])
    TorData=np.zeros((M,3))#Tcont,Rcol,Tor
    AB=['a','b']
    Tor=0
    Time=np.zeros((1),float)
    for i in range(0,M):
        [TorData[i,1], TorData[i,0]] = RcolwoATcontact(CrossData[2,i], CrossData[3,i], AB[int(CrossData[5,i])])
        TorData[i,2]=TorData[i,0]*TorData[i,1]*1e17 #Tor Time
        Tor=+TorData[i,2]
        StartTime=CrossData[0,i]-(TorData[i,0]/2)
        StopTime=CrossData[0,i]+(TorData[i,0]/2)
        numbertorsteps=np.floor(TorData[i,2]/torstep)

        Mat=np.linspace(StartTime,StopTime,numbertorsteps)
        Time = np.concatenate((Time, Mat))#,axis=1)
    '''
    print('no of torsteps',numbertorsteps)
    print(Time)
    print('Timelegnth',np.size(Time))
    print('TorDistLegnth',np.size(TorDist[0,0,:]))
    '''
    #plt.figure()
    plt.plot(Time,(TorDist[2,0,0:np.size(Time)]/(TorDist[0,0,0]+TorDist[1,0,0])),label='Precessing')
    plt.xlabel('Time yrs')
    plt.ylabel('Fractional Vol in Child')
    plt.legend()
    return


def analyticalevolutionMono(rfrag):
    lambdav=lambdatimeaverage(1e4)
    nini=TotalVol/(((4*pi)/3)*(rfrag**3.))
    N=int(1e4)
    Time=np.linspace(0,1e8,N)
    FracinChild=np.zeros(N)
    for i in range (0,N):
        FracinChild[i]=(4*nini*lambdav*Time[i])/(1+(4*nini*lambdav*Time[i]))

    plt.plot(Time,FracinChild,label='Mono Analytic')
    plt.legend()
    return





def analyticalevoltuionDistribution():
    K=(3/(8*pi))*(TotalVol/(rmax**0.5))
    nminin=K*(rmin**-3.5)
    nmaxin=K*(rmax**-3.5)
    lambdav = lambdatimeaverage(1e4)
    N=int(1e4)
    Time = np.linspace(0, 1e8, N)
    ParentFrac=np.zeros((2,N))
    for i in range(0,N):
        ParentFrac[0,i]=((1+(4*nminin*lambdav*Time[i])))**-(((rmax/rmin)**2.)/4) #parent
        ParentFrac[1,i]=1-ParentFrac[0,i] #child

    '''
    plt.figure()
    plt.plot(Time,ParentFrac,label='Analytical Evolution')
    plt.xlabel('Tor')
    plt.ylabel('scaled analytical Mass')
    '''

    #plt.figure()
    plt.plot(Time,ParentFrac[1,:] , label='Analytical Evolution')
    #plt.xlabel('Tor')
    #plt.ylabel('scaled analytical Child Mass')
    plt.legend()
    return


#PrecessingEvolutionGraphs()

#SaveTorDist(100)
#TorMassChange()
#PrecessingEvolutionGraphs()
#analyticalevolutionMono(1)
#analyticalevoltuionDistribution()
plt.show()
IniD=InitialDistCalc(5)
print(IniD)



