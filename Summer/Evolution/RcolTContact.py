from Functions import *
'''
au = 149597871e3
G = 6.67408e-11
Mstar = 1.2e30
'''

#(a1,e1)=orbitalvalues(3)
#(a2,e2)=orbitalvalues(10)
#a1=0.999*a1

#Values
'''
a1 = 2.5 * au
e1 = 0.998
a2 = 2.4 * au
e2 = 0.997
'''

#I=1* ((2 * pi) / 360)
#Rbeam=100e3

def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)*(10**6.) #in Yrs
    wp = (2 * np.pi) / Tp
    return (Tp,wp)




def RcolwoATcontact(a1,e1,a2,e2,s1,s2,x,I1,I2,Rbeam):
    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)
    CollisionData=CollisionPoints(a1,e1,s1,a2,e2,s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    if R==0:
        RColwoA=0
        Tcontact=0
        #print('no contact')

    else:
        td1=thetadot(a1,e1,C-s1)
        rd1=rdot(a1,e1,C-s1)
        V1=np.array([rd1,R*td1*cos(I1),R*td1*sin(I1)])
        v1=np.sqrt((rd1**2.)+((R*td1)**2.))
        td2 = thetadot(a2, e2, C - s2)
        rd2 = rdot(a2, e2, C - s2)
        V2 = np.array([rd2, R * td2 * cos(I2), R * td2 * sin(I2)])
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



        RColwoA=(16 / (3 * pi)) * (1 / (T1 * T2)) * ((vrelpar )/ (sinval*v1 * v2 ))*(year/Rbeam)  #

        Tcontact=abs((4*Rbeam*np.linalg.norm(np.cross(V1, V2)))/(((R** 2.) * td1 * td2 * abs(sin(I2-I1)))*(((wp1) * (rd1 / td1)) - (wp2) * (rd2 / td2))))


    return [RColwoA,Tcontact]



def RcolwoATcontactSecondOrder(a1,e1,a2,e2,s1,s2,x,I1,I2,Rbeam):
    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)
    CollisionData=CollisionPoints(a1,e1,s1,a2,e2,s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]
    else:
        print('Error, no point selected')
        R=0
    if R==0:
        #print('no  contact error')
        RColwoA=0
        Tcontact=0
    else:
        td1=thetadot(a1,e1,C-s1)
        rd1=rdot(a1,e1,C-s1)
        V1=np.array([rd1,R*td1*cos(I1),R*td1*sin(I1)])
        v1=np.sqrt((rd1**2.)+((R*td1)**2.))
        td2 = thetadot(a2, e2, C - s2)
        rd2 = rdot(a2, e2, C - s2)
        V2 = np.array([rd2, R * td2 * cos(I2), R * td2 * sin(I2)])
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



        RColwoA=(16 / (3 * pi)) * (1 / (T1 * T2)) * ((vrelpar )/ (sinval*v1 * v2 ))*(year/Rbeam)  #



        dr=((2*Rbeam)*np.linalg.norm(np.cross(V1, V2)))/((R** 2.) * td1 * td2 * abs(sin(I2-I1)))
        Alpha=(wp1 * (rd1 / td1)) - (wp2 * (rd2 / td2))
        Beta=abs(((((wp1 * (rd1 / td1))**2.) - ((wp2 * (rd2 / td2))**2.))/R)+((R/2)*((wp1**2.)-(wp2**2.)))-(((R**2.)/2)*((((wp1**2.)/(a1*(1-(e1**2.))))-((wp2**2.)/(a2*(1-(e2**2.))))))))
        Tcontact = (np.sqrt((Alpha ** 2.) + (4 * Beta * dr)) - np.sqrt((Alpha ** 2.) - (4 * Beta * dr))) / (2*Beta)


    return [RColwoA,Tcontact]


def TcontactSearchMethod(a1,e1,a2,e2,s1i,s2i,x,I1,I2,Rbeam):

    CollisionData=CollisionPoints(a1,e1,s1i,a2,e2,s2i)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]
    else:
        print('Error, no point selected')
        R=0
    if R==0:
        #print('no  contact error')
        RColwoA=0
        Tcontact=0
        dr=0

    else:
        td1=thetadot(a1,e1,C-s1i)
        rd1=rdot(a1,e1,C-s1i)
        V1=np.array([rd1,R*td1*cos(I1),R*td1*sin(I1)])
        td2 = thetadot(a2, e2, C - s2i)
        rd2 = rdot(a2, e2, C - s2i)
        V2 = np.array([rd2, R * td2 * cos(I2), R * td2 * sin(I2)])
        dr = ((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1)))
        tinterval=0.5

        (Tp1, wp1) = TPrecess(a1, e1)
        (Tp2, wp2) = TPrecess(a2, e2)
        # up loop
        dr = [((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1))),0]
        dRplus=[0,0]
        y=0
        tplus=0
        while abs(dRplus[y])<dr[y]:
            tplus += tinterval
            y=(y+1)%2
            dRplus[y] = npr(a1, e1, C - s1i - (wp1 * tplus)) - npr(a2, e2, C - s2i - (
        wp2 * tplus))  # changed from plue to minus!
            td1 = thetadot(a1, e1, C - s1i - (wp1 * tplus))
            rd1 = rdot(a1, e1, C - s1i - (wp1 * tplus))
            V1 = np.array([rd1, R * td1 * cos(I1), R * td1 * sin(I1)])
            td2 = thetadot(a2, e2, C - s2i - (
        wp2 * tplus))
            rd2 = rdot(a2, e2, C - s2i - (
        wp2 * tplus))
            V2 = np.array([rd2, R * td2 * cos(I2), R * td2 * sin(I2)])

            dr[y] = ((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1)))
        if abs(abs(dRplus[y])-dr[y])>abs(abs(dRplus[(y+1)%2])-dr[(y+1)%2]): #if the previous was closer, go back
            tplus-=tinterval

        dRminus = [0, 0]
        dr = [((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1))), 0]
        y = 0
        tminus = 0
        while abs(dRminus[y]) < dr[y]:
            tminus -= tinterval
            y = (y + 1) % 2
            dRminus[y] = npr(a1, e1, C - s1i - (wp1 * tminus)) - npr(a2, e2, C - s2i - (
                wp2 * tminus))
            td1 = thetadot(a1, e1, C - s1i - (wp1 * tminus))
            rd1 = rdot(a1, e1, C - s1i - (wp1 * tminus))
            V1 = np.array([rd1, R * td1 * cos(I1), R * td1 * sin(I1)])
            td2 = thetadot(a2, e2, C - s2i - (
                wp2 * tminus))
            rd2 = rdot(a2, e2, C - s2i - (
                wp2 * tminus))
            V2 = np.array([rd2, R * td2 * cos(I2), R * td2 * sin(I2)])

            dr[y] = ((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1)))
        if abs(abs(dRminus[y]) - dr[y]) > abs(
                        abs(dRminus[(y + 1) % 2]) - dr[(y + 1) % 2]):  # if the previous was closer, go back
            tminus += tinterval
        '''
        dRminus = 0
        tminus = 0
        dr = ((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1)))
        while abs(dRminus) < dr:
            tminus -= tinterval
            dRminus = npr(a1, e1, C- s1i - (wp1 * tminus)) - npr(a2, e2, C - s2i - (
                wp2 * tminus))  # changed from plue to minus!
            td1 = thetadot(a1, e1, C - s1i - (wp1 * tplus))
            rd1 = rdot(a1, e1, C - s1i - (wp1 * tplus))
            V1 = np.array([rd1, R * td1 * cos(I1), R * td1 * sin(I1)])
            td2 = thetadot(a2, e2, C - s2i - (
                wp2 * tplus))
            rd2 = rdot(a2, e2, C - s2i - (
                wp2 * tplus))
            V2 = np.array([rd2, R * td2 * cos(I2), R * td2 * sin(I2)])

            dr = ((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1)))
        '''


        #print(tplus)
        #print(tminus)
        Tcontact=abs(tplus-tminus)
        dr = ((2 * Rbeam) * np.linalg.norm(np.cross(V1, V2))) / ((R ** 2.) * td1 * td2 * abs(sin(I2 - I1)))


    return [Tcontact,dr]


def TorGraph():
    N=int(1000)
    s1=0
    s2=np.linspace(0,2*pi,N)
    Data=np.zeros((3,N))#Rcol,Tcont,tor
    for i in range(1,N):
        [Data[0,i],Data[1,i]]=RcolwoATcontact(0,s2[i],'b')
        Data[2,i]=Data[0,i]*Data[1,i]*1e17

    plt.figure()
    plt.semilogy(s2[1:N-1],Data[0,1:N-1],label='Rcol')
    plt.xlabel('Lambda')
    plt.ylabel('Rcol')

    plt.figure()
    plt.semilogy(s2[1:N-1], Data[1, 1:N-1], label='Tcontact')
    plt.xlabel('Lambda')
    plt.ylabel('Tcontact, yrs')

    plt.figure()
    plt.semilogy(s2[1:N-1], Data[2, 1:N-1], label='Tor x 10^17')
    plt.xlabel('Lambda')
    plt.ylabel('Tor')
    return

def lambdatimeaverage(N):
    N=int(N)
    #N = int(1000)
    s1 = 0
    s2 = np.linspace(0, 2 * pi, N)
    Data = np.zeros((3, N))  # Rcol,Tcont,tor
    LambAv=0
    for i in range(1, N):
        [Data[0, i], Data[1, i]] = RcolwoATcontact(0, s2[i], 'b')
        Data[2, i] = Data[0, i] * Data[1, i]
        LambAv +=Data[2, i]
        [Data[0, i], Data[1, i]] = RcolwoATcontact(0, s2[i], 'a')
        Data[2, i] = Data[0, i] * Data[1, i]
        LambAv+=Data[2,i]

    Tpmin=min(Tp1,Tp2)
    LambAv = (LambAv / (2 * N))*(4/Tpmin)
    print(LambAv)

    return LambAv

def TcontactGraphs(N,rs1,rs2,pmod):
    N=int(N)
    Rbeam=100*1e3
    I1=0
    I2=10*(pi/180)

    (a1, e1) = orbitalvalues(rs1)
    (a2, e2) = orbitalvaluesmod(rs2, pmod)



    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)

    L = np.linspace(0, 2 * pi, N)  # NOTE CHECK L def!! + or - #s1-s2
    R = np.zeros((2, N))  # radial collision points
    C = np.zeros((2, N))  # theta collision point
    AB = ['a', 'b']
    # s

    for i in range(1, N):
        [[C[0, i], C[1, i]], [R[0, i], R[1, i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    Rdot1 = rdot(a1, e1, C[:, :])
    Thetadot1 = thetadot(a1, e1, C[:, :])
    Rdot2 = rdot(a2, e2, C[:, :] - L)
    Thetadot2 = thetadot(a2, e2, C[:, :] - L)
    Vinc1 = np.zeros((2, 3, N))
    Vinc2 = np.zeros((2, 3, N))
    vrelpar = np.zeros((2, N))



    Tcontact= np.zeros((2, N))
    TcontactSecond = np.zeros((2, N))
    TcontSearch= np.zeros((2, N))
    Rcol= np.zeros((2, N))
    ProblemValue=np.zeros((2, N))
    dr=np.zeros((2, N))
    print('Calculating times...')
    for i in range(1, N - 1):
        #print(i,'/',N-2)
        for x in (0, 1):
            #MinTpr[x,i] = TPrAnalytic(OrbitdataVdata[x, 0, i], OrbitdataVdata[x, 1, i])  #Min

            #MaxTpr[x,i] = MinTpr[x,i] * dflr[x, i]*rfrag  #Max
            [Rcol[x,i],Tcontact[x,i]]=RcolwoATcontact(a1,e1,a2,e2,0,L[i],AB[x],I1,I2,Rbeam)


            TcontactSecond[x, i] = RcolwoATcontactSecondOrder(a1, e1, a2, e2, 0, L[i], AB[x], I1, I2, Rbeam)[1]
            ProblemValue[x,i]=((Rdot1[x,i]/Thetadot1[x,i])*wp1)-((Rdot2[x,i]/Thetadot2[x,i])*wp2)
            [TcontSearch[x,i],dr[x,i]]=TcontactSearchMethod(a1, e1, a2, e2, 0, L[i], AB[x], I1, I2, Rbeam)

    print('plotting')
    for x in (0, 1):
        plt.figure()
        plt.title('Times at %s pts, rsource1=%s au, rsource2=%s au with pmod=%s' % (AB[x], rs1, rs2, pmod))
        plt.semilogy(L[1:N - 1] / pi, Tcontact[x, 1:N - 1], label='Contact time')
        plt.semilogy(L[1:N - 1] / pi, TcontactSecond[x, 1:N - 1], label='Second Order Contact time')
        plt.semilogy(L[1:N - 1] / pi, TcontSearch[x, 1:N - 1], label='Searched Tcont')
        plt.ylabel('Times, yrs')
        plt.xlabel('Lambda')
        plt.legend()
    '''
    plt.figure()
    plt.title('dr at %s pts, rsource1=%s au, rsource2=%s au with pmod=%s' % (AB[x], rs1, rs2, pmod))
    plt.semilogy(L[1:N - 1] / pi, dr[x, 1:N - 1], label='dr')
    plt.ylabel('dr, distance')
    plt.xlabel('Lambda')
    plt.legend()

    plt.figure()
    plt.title('ProbValue at %s pts, rsource1=%s au, rsource2=%s au with pmod=%s' % (AB[x], rs1, rs2, pmod))
    plt.semilogy(L[1:N - 1] / pi, abs(ProblemValue[x, 1:N - 1]), label='Prob Value')
    '''

    return




#COMMANDS
#[RColwoA,Tcontact]=RcolwoATcontact(0,3,'a')
#print([RColwoA,Tcontact])
#TorGraph()
#lambdatimeaverage(1e4)

TcontactGraphs(1e3,10,3,0.9)
plt.show(), print('plotted')

#3,10,0.9