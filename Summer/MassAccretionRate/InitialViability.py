#from Evolution.RcolTContact import *
from Functions import *

def MonoGraphs(): #Plots initial rates for fragment and parent vols
    #Target Rates in gs per s
    N=100
    Mdot=np.zeros((6))
    rsource1 = 3
    rsource2 = 3.1
    pmod = 0.9
    (a1, e1) = orbitalvalues(rsource1)
    (a2, e2) = orbitalvaluesmod(rsource2, pmod)


    Kappa=lambdatimeaverage(a1,e1,a2,e2,int(1e4))*(6/pi)*density*1e3  #changed to grams!
    ''''
    V=np.linspace(12,18,N)
    V=np.power(10,V)
    r=np.zeros((6,N))
    plt.figure()
    for i in range(0,6):
        Mdot[i]=10**(i+6)
        for j in range(0,N):
            r[i,j]=(Kappa*(V[j]**2.))/(Mdot[i]*year)


        plt.loglog(V,r[i,:],label='Mdot=10^%s'%(i+6))

    plt.ylabel('frag radius, m')
    plt.xlabel('Parent Volumes, m^3')
    plt.title(' rs1=%s au, rs2=%s au with pmod=%s'%(rsource1,rsource2,pmod))
    #plt.ylim([10**-2,10**2])
    plt.legend()
    '''

    Rast = np.linspace(3, 6, N)
    Rast = np.power(10, Rast)
    r = np.zeros((6, N))
    plt.figure()
    for i in range(0, 6):
        Mdot[i] = 10 ** (i + 6)
        for j in range(0, N):
            r[i, j] = (Kappa * (Rast[j] ** 6.)*(((4/3)*pi))**2.) / (Mdot[i] * year)

        plt.loglog(Rast, r[i, :], label='Mdot=10^%s' % (i + 6))

    plt.ylabel('frag radius, m')
    plt.xlabel('Rast, m')
    plt.title(' rs1=%s au with p1=%s au, rs2=%s au with p2=%s au' % (rsource1,0.01, rsource2, 0.009))#pmod*0.01))
    # plt.ylim([10**-2,10**2])
    plt.legend()


    return


I=1* ((2 * pi) / 360)
Rbeam=100e3

def TPrecess(a,e):
    #print(a/au,e)
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)*(10**6.) #in Yrs
    #print(Tp)
    wp = (2 * np.pi) / Tp
    return (Tp,wp)




def RcolwoATcontactParam(a1,e1,wp1,s1,a2,e2,wp2,s2,x):
    CollisionData=CollisionPoints(a1,e1,s1,a2,e2,s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]
    else:
        print('Error, no point selected')
    if R==0:
        RColwoA=0
        Tcontact=0

    else:
        td1=thetadot(a1,e1,C-s1)
        rd1=rdot(a1,e1,C-s1)
        V1=np.array([rd1,R*td1*cos(I),R*td1*sin(I)])
        v1=np.sqrt((rd1**2.)+((R*td1)**2.))
        td2 = thetadot(a2, e2, C - s2)
        rd2 = rdot(a2, e2, C - s2)
        V2 = np.array([rd2, R * td2 * cos(I), R * td2 * sin(I)])
        v2 = np.sqrt((rd2 ** 2.) + ((R * td2) ** 2.))
        Vrelpar=V1-V2
        vrelpar=npal.norm(Vrelpar)
        gammangle=np.arccos(np.vdot(V1,V2)/(v1*v2))


        T1 = 2 * pi * np.sqrt((a1 ** 3.) / (G * Mstar))  # Orbital Period
        T2 = 2 * pi * np.sqrt((a2 ** 3.) / (G * Mstar))

        sinval=abs(sin(gammangle))
        '''
        if (16*abs((rd1/td1)-(rd2/td2)))/(3*pi*R)>sinval:
            sinval=(16*abs((rd1/td1)-(rd2/td2)))/(3*pi*R)
            print('sinval',sinval)
            print('lambda=',(s1-s2)%(2*pi))
        '''



        RColwoA=(16 / (3 * pi)) * (1 / (T1 * T2)) * ((vrelpar )/ (sinval*v1 * v2 ))*(year/Rbeam)*0.46  #

        Tcontact=abs((4*Rbeam*np.linalg.norm(np.cross(V1, V2)))/(((R** 2.) * td1 * td2 * sin(I))*((wp1) * (rd1 / td1)) - (wp2) * (rd2 / td2)))


    return [RColwoA,Tcontact]
def lambdatimeaverage(a1,e1,a2,e2,N):
    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)
    #N=int(N)
    #N = int(1000)
    s1 = 0
    s2 = np.linspace(0, 2 * pi, N)
    Data = np.zeros((3, N))  # Rcol,Tcont,tor
    LambAv=0
    for i in range(1, N):
        [Data[0, i], Data[1, i]] = RcolwoATcontactParam(a1,e1,wp1,0,a2,e2,wp2, s2[i], 'b')
        Data[2, i] = Data[0, i] * Data[1, i]
        LambAv +=Data[2, i]
        [Data[0, i], Data[1, i]] = RcolwoATcontactParam(a1,e1,wp1,0,a2,e2,wp2, s2[i], 'a')
        Data[2, i] = Data[0, i] * Data[1, i]
        LambAv+=Data[2,i]

    Tsep=min(Tp1,Tp2)/4
    LambAv = (LambAv / (2 * (N-2)))/Tsep
    #print(LambAv)

    return LambAv

def paramterlamdagraph():
    N = 200
    rsource = np.linspace(2, 15, N)
    Orbits1 = np.zeros((2,N))
    Orbits2 = np.zeros((2, N))
    for i in range(0,N):
        Orbits1[:,i]=orbitalvalues(rsource[i])
        Orbits2[:, i] = orbitalvaluesmod(rsource[i],0.9)
    #Orbits[0, :] = Orbits[0, :] * au
    LambData = np.zeros((N, N))
    Rast=1e5
    Rfrag=1
    Nini=(Rast/Rfrag)**3.
    Thalf = np.zeros((N, N))

    for i in range(0, N):
        print(i,'/%s'%N)
        for j in range(0, N):
            #print(Orbits[0,i]/au,Orbits[1,i])
            LambData[i,j]=np.log10(lambdatimeaverage(Orbits1[0,i],Orbits1[1,i],Orbits2[0,j],Orbits2[1,j],360))
            Thalf[i,j]=-LambData[i,j]-np.log10(Nini)

    print('LambavFound')
    plt.figure()

    plt.title('log10(xi)')
    cp = plt.contour(rsource, rsource, LambData)#, np.arange(0, 1.2, 0.2))  # , colors='k')
    plt.clabel(cp, inline=1, fontsize=10)
    plt.xlabel('rsource au, pmod0.9')
    plt.ylabel('rsource au')

    plt.figure()
    plt.title('log10(Thalf/years) ')
    cp = plt.contour(rsource, rsource, Thalf)  # , np.arange(0, 1.2, 0.2))  # , colors='k')
    plt.clabel(cp, inline=1, fontsize=10)
    plt.xlabel('rsource au, pmod0.9')
    plt.ylabel('rsource au')
    return


def MonoLambavRastGraphs():
    N=1000
    ttop=11
    tmin=5
    rfrag=1
    tlegnth=ttop-tmin
    iMdot = np.zeros((tlegnth))  # g/s
    lambdav=np.zeros((tlegnth,N))

    Rast=np.linspace(4,7,N)
    plt.figure()
    for i in range(0, N):
        Rast[i]=10**Rast[i]
    for t in range(0,ttop-tmin):
        iMdot[t]=10**(t+tmin)
        for i in range(0,N):
            lambdav[t,i]=(9/(64*(pi**2.)*density))*(rfrag/(Rast[i]**6.))*iMdot[t]*(year/1000)

        plt.loglog(Rast,lambdav[t,:],label='Mdot=10^%s g/s'%(t+tmin))
    plt.legend()
    plt.xlabel('Rast, m')
    plt.ylabel('xi')

    lambdav=np.linspace(20,23,N)
    Rast=np.zeros((tlegnth,N))
    plt.figure()
    for i in range(0, N):
        lambdav[i] = 10 ** (-lambdav[i])
    for t in range(0, ttop - tmin):
        iMdot[t] = 10 ** (t + tmin)
        for i in range(0, N):
            Rast[t, i] = ((9 / (64 * (pi ** 2.) * density)) * (rfrag / (lambdav[i] )) * iMdot[t] * (year / 1000))**(1/6)

        plt.loglog(Rast[t, :], lambdav, label='Mdot=10^%s g/s' % (t + tmin))
    plt.legend()
    plt.xlabel('Rast, m')
    plt.ylabel('xi')
    plt.title('Graph of Mdot dependancies on Rast and xi')

    return
#paramterlamdagraph()
#MonoLambavRastGraphs()
MonoGraphs()
plt.show()




