from Functions import *
from Summer.PR.TprAnalytic import *
from Summer.Evolution.RcolTContact import *
rsource1=3
rsource2=10
pmod=0.9
s1i=1
s2i=0
(ap1,ep1)=orbitalvalues(rsource1)
(ap2,ep2)=orbitalvaluesmod(rsource2,pmod)
I1=0
I2 = 1 * ((2 * pi) / 360)
rCollision  =0.1*au
Rbeam=100*1e3

def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)*(10**6.) #in Yrs
    wp = (2 * np.pi) / Tp
    return (Tp,wp)

def RcolwoATcontactVelocities(rd1,td1,a1,rd2,td2,a2,R,wp1,wp2, Inc1,Inc2):


        V1 = np.array([rd1, R * td1 * cos(Inc1), R * td1 * sin(Inc1)])
        v1 = np.sqrt((rd1 ** 2.) + ((R * td1) ** 2.))
        V2 = np.array([rd2, R * td2 * cos(Inc2), R * td2 * sin(Inc2)])
        v2 = np.sqrt((rd2 ** 2.) + ((R * td2) ** 2.))
        Vrelpar = V1 - V2
        vrelpar = npal.norm(Vrelpar)
        gammangle = np.arccos(np.vdot(V1, V2) / (v1 * v2))

        T1 = 2 * pi * np.sqrt((a1 ** 3.) / (G * Mstar))  # Orbital Period
        T2 = 2 * pi * np.sqrt((a2 ** 3.) / (G * Mstar))

        sinval = abs(sin(gammangle))
        '''
        if (16*abs((rd1/td1)-(rd2/td2)))/(3*pi*R)>sinval:
            sinval=(16*abs((rd1/td1)-(rd2/td2)))/(3*pi*R)
            print('sinval',sinval)
            print('lambda=',(s1-s2)%(2*pi))
        '''

        RColwoA = (16 / (3 * pi)) * (1 / (T1 * T2)) * ((vrelpar) / (sinval * v1 * v2)) * (year / Rbeam)  #

        Tcontact = abs((4 * Rbeam * np.linalg.norm(np.cross(V1, V2))) / (
        ((R ** 2.) * td1 * td2 * abs(sin(Inc1-Inc2))) * ((wp1) * (rd1 / td1)) - (wp2) * (rd2 / td2)))

        return [RColwoA, Tcontact]


def PRVsp():
    print('Creating Data...')
    N=int(1e3)
    rd=np.linspace(-150e3,150e3,2*N)
    tdr=np.linspace(1,100e3,N)
    semimajora=np.zeros((2*N,N))
    eccentricity=np.zeros((2*N,N))
    Tpr = np.zeros((2 * N, N))
    for i in range(0,N):

        for j in range(0,2*N):

            eccentricity[j, i] = newe(rd[j], tdr[i] / rCollision, rCollision)
            if eccentricity[j,i]<0.99:
                semimajora[j, i] = newa(rd[j], tdr[i] / rCollision, rCollision)
                Tpr[j,i]=log10(TPrAnalytic(semimajora[j,i],eccentricity[j,i]))
    print('Plotting...')

    #Contour
    plt.figure()
    plt.title('Eccentricity')
    cp=plt.contour(rd,tdr,np.transpose(eccentricity),np.arange(0,1.2,0.2))#, colors='k')
    plt.clabel( cp,inline=1, fontsize=10)
    plt.figure()
    plt.title('SemiMajor Axis/au')
    semia = plt.contour(rd, tdr, np.transpose(semimajora/au),np.arange(0,5,0.2))#, colors='k')
    plt.clabel(semia, inline=1, fontsize=10)
    plt.figure()
    plt.title('Tpr, yrs')
    Tprgraph = plt.contour(rd, tdr, np.transpose(Tpr),np.arange(0,10,1), label='PRTimescale')
    plt.clabel(Tprgraph, inline=1, fontsize=10)
    cs2 = plt.contour(Tprgraph, levels=Tprgraph.levels[::2],colors='r')
    '''

    #3d
    eccengraph= plt.figure()
    plt.title('Eccentricity')
    ax = eccengraph.add_subplot(111, projection='3d')
    ax.plot_surface(rd, tdr, np.transpose(eccentricity), color='b')
    '''
    print('Plotted')
    return


def VspHValues(s1, s2, x):
    N = 1000

    CollisionData = CollisionPoints(ap1, ep1, s1, ap2, ep2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]
    else:
        print('Error, no point selected')
    td1 = thetadot(ap1, ep1, C - s1)
    rd1 = rdot(ap1, ep1, C - s1)
    V1 =  np.array([rd1, R * td1 * cos(I1), R * td1 * sin(I1)])#np.array([rd1, R * td1])  #
    v1 = np.sqrt((rd1 ** 2.) + ((R * td1) ** 2.))
    td2 = thetadot(ap2, ep2, C - s2)
    rd2 = rdot(ap2, ep2, C - s2)
    V2 =  np.array([rd2, R * td2 * cos(I2), R * td2 * sin(I2)])# np.array([rd2, R * td2])  #
    v2 = np.sqrt((rd2 ** 2.) + ((R * td2) ** 2.))
    Vrelpar = V2-V1
    vrelpar = npal.norm(Vrelpar)
    gammangle = np.arccos(np.vdot(V1, V2) / (v1 * v2))

    Velocities = np.zeros((3, N))
    H = np.linspace(0, 1, N)


    Data = np.zeros((11,N)) #a,e,i,Tpr,appo,Tprecess,wp,Tcont1,Tcont2,Rcol1,Rcol2
    Vdata = np.zeros(( 2, N))  # rd,Rtd

    for i in range(0, N):
        Velocities[:, i] = V1 + (H[i] * Vrelpar)
        Vdata[ :, i] = np.array([Velocities[ 0, i], np.sqrt((Velocities[ 1, i] ** 2.) + (Velocities[ 2, i] ** 2.))]) #rd,R*td
    for i in range(0,N):
        Data[0,i]=newa(Vdata[0,i],Vdata[1,i]/R,R) #a
        Data[1, i] = newe(Vdata[0, i], Vdata[1, i] / R,R) #e
        Data[2,i]= np.arccos(Velocities[ 1, i] / Vdata[ 1, i]) #Inclination
        if Data[2,i]>pi:
            Data[2,i]=(2*pi)-Data[2,i]
        Data[3,i] = TPrAnalytic(Data[0,i],Data[1,i]) #Tpr
        Data[4,i]=Data[0,i]*(1+Data[1,i]) # appocentre
        [Data[5,i],Data[6, i]]=TPrecess(Data[0,i],Data[1,i]) #Tprecess,wp
    Data[2,0]=I1

    for i in range(1, N):
        [Data[9, i],Data[7,i]] = RcolwoATcontactVelocities(Vdata[0, 0], Vdata[1, 0] / R,Data[0,0], Vdata[0, i], Vdata[1, i] / R,Data[0,i], R, Data[6, 0],
                                  Data[6, i],Data[2,0],Data[2,i])  # Rcol1, Tcont1

    for i in range(0, N-1):
        [Data[10, i],Data[8, i]] = RcolwoATcontactVelocities(Vdata[0,i],Vdata[1,i]/R,Data[0,i],Vdata[0,N-1],Vdata[1,N-1]/R,Data[0,N-1],R,Data[6,i],Data[6, N-1],Data[2,i],Data[2,N-1]) #Tcont2

    #Plots
    #a Plot
    plt.figure()
    plt.semilogy(H,Data[0,:]/au,label='Semimajor Axis')
    plt.xlabel('H')
    plt.ylabel('Semimajor Axis au')

    #e label
    plt.figure()
    plt.plot(H, Data[1, :], label='Eccentricity')
    plt.xlabel('H')
    plt.ylabel('Eccentricity')


    #Appocentre label
    plt.figure()
    plt.semilogy(H,Data[4,:]/au , label='Appocentre')
    plt.semilogy(np.array([0,1]), np.array([R,R])/au, label='Collision Pt')
    plt.xlabel('H')
    plt.ylabel('Appocentre')

    #Tpr
    plt.figure()
    plt.semilogy(H, Data[3, :], label='Tpr')
    plt.xlabel('H')
    plt.ylabel('Tpr, years')
    plt.title('Tpr against H at rs1=%s au, rs2=%s au with pmod=%s, Lambda=%s'%(rsource1,rsource2,pmod,s2-s1))
    # TPrecess
    plt.figure()
    plt.semilogy(H, Data[5, :], label='TPrecess')
    plt.xlabel('H')
    plt.ylabel('TPrecess, years')
    plt.title('TPrecess against H at rs1=%s au, rs2=%s au with pmod=%s, Lambda=%s' % (rsource1, rsource2, pmod, s2 - s1))

    # Tcontact
    plt.figure()
    plt.semilogy(H, Data[7, :], label='TCont P1')
    plt.semilogy(H, Data[8, :], label='TCont P2')
    plt.xlabel('H')
    plt.ylabel('TCont, years')
    plt.title('TCont against H at rs1=%s au, rs2=%s au with pmod=%s, Lambda=%s' % (rsource1, rsource2, pmod, s2 - s1))
    plt.legend()

    # Rcol
    plt.figure()
    plt.semilogy(H, Data[9, :], label='Rcol 1')
    plt.semilogy(H, Data[10, :], label='Rcol 2')
    plt.xlabel('H')
    plt.ylabel('Rcol')
    plt.title('Rcol against H at rs1=%s au, rs2=%s au with pmod=%s, Lambda=%s' % (rsource1, rsource2, pmod, s2 - s1))
    plt.legend()



    return
VspHValues(0,2,'a')

#PRVsp()

plt.show()
