from Functions import *
from Summer.PR.TprAnalytic import *
rsource1=3
rsource2=10
pmod=0.9
s1i=1
s2i=0
(a1,e1)=orbitalvalues(rsource1)
(a2,e2)=orbitalvaluesmod(rsource2,pmod)
I = 1 * ((2 * pi) / 360)
R=0.1*au

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

            eccentricity[j, i] = newe(rd[j], tdr[i] / R, R)
            if eccentricity[j,i]<0.99:
                semimajora[j, i] = newa(rd[j], tdr[i] / R, R)
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

    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]
    else:
        print('Error, no point selected')
    td1 = thetadot(a1, e1, C - s1)
    rd1 = rdot(a1, e1, C - s1)
    V1 = np.array([rd1, R * td1])  # np.array([rd1, R * td1 * cos(I), R * td1 * sin(I)])
    v1 = np.sqrt((rd1 ** 2.) + ((R * td1) ** 2.))
    td2 = thetadot(a2, e2, C - s2)
    rd2 = rdot(a2, e2, C - s2)
    V2 = np.array([rd2, R * td2])  # np.array([rd2, R * td2 * cos(I), R * td2 * sin(I)])
    v2 = np.sqrt((rd2 ** 2.) + ((R * td2) ** 2.))
    Vrelpar = V2-V1
    vrelpar = npal.norm(Vrelpar)
    gammangle = np.arccos(np.vdot(V1, V2) / (v1 * v2))

    Momentas = np.zeros((2, N))
    H = np.linspace(0, 1, N)

    Data = np.zeros((4,N)) #a,e,Tpr

    for i in range(0, N):
        Momentas[:, i] = V1 + H[i] * Vrelpar
        Data[0,i]=newa(Momentas[0,i],Momentas[1,i]/R,R) #a
        Data[1, i] = newe(Momentas[0, i], Momentas[1, i] / R,R) #e
        Data[2,i] = TPrAnalytic(Data[0,i],Data[1,i]) #Tpr
        Data[3,i]=Data[0,i]*(1+Data[1,i]) # appocentre
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
    plt.semilogy(H,Data[3,:]/au , label='Appocentre')
    plt.semilogy(np.array([0,1]), np.array([R,R])/au, label='Collision Pt')
    plt.xlabel('H')
    plt.ylabel('Appocentre')

    #Tpr
    plt.figure()
    plt.semilogy(H, Data[2, :], label='Tpr')
    plt.xlabel('H')
    plt.ylabel('Tpr, years')


    return
VspHValues(0,1,'a')

#PRVsp()

plt.show()
