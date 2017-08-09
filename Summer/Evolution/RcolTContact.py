from Functions import *
'''
au = 149597871e3
G = 6.67408e-11
Mstar = 1.2e30
'''

(a1,e1)=orbitalvalues(3)
(a2,e2)=orbitalvalues(10)
a1=0.999*a1

#Values
'''
a1 = 2.5 * au
e1 = 0.998
a2 = 2.4 * au
e2 = 0.997
'''

I=1* ((2 * pi) / 360)
Rbeam=100e3  #CHANGED

def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)*(10**6.) #in Yrs
    wp = (2 * np.pi) / Tp
    return (Tp,wp)

(Tp1,wp1)=TPrecess(a1,e1)
(Tp2,wp2)=TPrecess(a2,e2)


def RcolwoATcontact(a1,e1,a2,e2,s1,s2,x):
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



        RColwoA=(16 / (3 * pi)) * (1 / (T1 * T2)) * ((vrelpar )/ (sinval*v1 * v2 ))*(year/Rbeam)  #

        Tcontact=abs((4*Rbeam*np.linalg.norm(np.cross(V1, V2)))/(((R** 2.) * td1 * td2 * sin(I))*((wp1) * (rd1 / td1)) - (wp2) * (rd2 / td2)))


    return [RColwoA,Tcontact]

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



#COMMANDS
#[RColwoA,Tcontact]=RcolwoATcontact(0,3,'a')
#print([RColwoA,Tcontact])
#TorGraph()
#lambdatimeaverage(1e4)


plt.show()