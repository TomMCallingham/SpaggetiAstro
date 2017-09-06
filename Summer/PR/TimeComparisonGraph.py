from Functions import *
from PR.TprAnalytic import *
from Fragmentation.FragFuncs import *
from GRPrecession.Precession import *
from Evolution.RcolTContact import *
density=2000
Qa=620
Qb=(5e-6)*density
a=0.3
b=1.5
rsource1=10
rsource2=10
pmod=0.9
(a1,e1)=orbitalvalues(rsource1)
(a2,e2)=orbitalvaluesmod(rsource2,pmod)
#a1=0.999*a1
I2=1* ((2 * pi) / 360)
I1=0
Rbeam=100*1e3

def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)*(10**6.) #in Yrs
    wp = (2 * np.pi) / Tp
    return (Tp,wp)

(Tp1,wp1)=TPrecess(a1,e1)
(Tp2,wp2)=TPrecess(a2,e2)


def TimeGraphs(rfrag, N):
    # N = 1000  # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)  # NOTE CHECK L def!! + or - #s1-s2
    R = np.zeros((2, N))  # radial collision points
    C = np.zeros((2, N))  # theta collision point
    AB=['a','b']
    # s
    print('Creating Parent Values...')
    for i in range(1, N):
        [[C[0, i], C[1, i]], [R[0, i], R[1, i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    Rdot1 = rdot(a1, e1, C[:, :])
    Thetadot1 = thetadot(a1, e1, C[:, :])
    Rdot2 = rdot(a2, e2, C[:, :] - L)
    Thetadot2 = thetadot(a2, e2, C[:, :] - L)
    Vinc1 = np.zeros((2, 3, N))
    Vinc2 = np.zeros((2, 3, N))
    vrelpar = np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            Vinc1[x, :, i] = np.array(
                [Rdot1[x, i], R[x, i] * Thetadot1[x, i] * cos(I1), R[x, i] * Thetadot1[x, i] * sin(I1)])
            Vinc2[x, :, i] = np.array(
                [Rdot2[x, i], R[x, i] * Thetadot2[x, i] * cos(I2), R[x, i] * Thetadot2[x, i] * sin(I2)])
            vrelpar[x, i] = np.linalg.norm(Vinc1[x, :, i] - Vinc2[x, :, i])
    print('Fragmenting...')
    # Fragmenting
    dflr = np.zeros((2, N))

    for i in range(1, N):
        dflr[0, i] = flrfunc(vrelpar[0, i], 2 * rfrag, 2 * rfrag) ** (1 / 3)
        dflr[1, i] = flrfunc(vrelpar[1, i], 2 * rfrag, 2 * rfrag) ** (1 / 3)
    # Find Collided Orbit Values
    print('Finding Child Orbits...')
    Vmid = (Vinc1 + Vinc2) / 2
    Vdata = np.zeros((2, 2, N))  # rd,Rtd
    OrbitdataVdata = np.zeros((2, 3, N))  # (a,e,i)
    for i in range(1, N - 1):
        for x in (0, 1):
            Vdata[x, :, i] = np.array([Vmid[x, 0, i], np.sqrt((Vmid[x, 1, i] ** 2.) + (Vmid[x, 2, i] ** 2.))])
            OrbitdataVdata[x, :, i] = np.array([newa(Vdata[x, 0, i], Vdata[x, 1, i] / R[x, i], R[x, i]),
                                                newe(Vdata[x, 0, i], Vdata[x, 1, i] / R[x, i], R[x, i]),
                                                np.arccos(Vmid[x, 1, i] / Vdata[x, 1, i])])
    print('Finding Times...')
    Tpr = np.zeros((2, N))
    Tcontact= np.zeros((2, N))
    for i in range(1, N - 1):
        for x in (0, 1):
            Tpr[x,i] = TPrAnalytic(OrbitdataVdata[x, 0, i], OrbitdataVdata[x, 1, i]) * dflr[x, i]
            Tcontact[x,i]=RcolwoATcontact(a1,e1,a2,e2,0,L[i],AB[x],I1,I2,Rbeam)[1]

    Tsep = min(Tp1,Tp2)/4

    for x in (0, 1):
        plt.figure()
        plt.title('Times at %s pts, rsource1=%s au, rsource2=%s au with pmod=%s'%(AB[x],rsource1,rsource2,pmod))
        plt.semilogy(L[1:N-1]/pi,Tpr[x,1:N-1],label='PR Time')
        plt.semilogy(L[1:N - 1]/pi, Tcontact[x, 1:N - 1], label='Contact time')
        plt.semilogy(np.array([0,2]), np.array([Tsep,Tsep]), label='Seperation Time')
        plt.ylabel('Times, yrs')
        plt.xlabel('Lambda')
        plt.legend()


    print('Times Found')

    # vrel plot
    plt.figure()
    plt.semilogy(L[1:N - 1], vrelpar[0, 1:N - 1], label='Vrel at a')
    plt.semilogy(L[1:N - 1], vrelpar[1, 1:N - 1], label='Vrel at b')
    plt.legend()
    plt.figure()
    plt.semilogy(L[1:N - 1], dflr[0, 1:N - 1], label='Reduction in Fragments')
    plt.semilogy(L[1:N - 1], dflr[1, 1:N - 1], label='Reduction in Fragments')
    plt.legend()



    return


TimeGraphs(1, 1000)
plt.show()

'''
colrad=0.2*au
rd=152949
trd=2228218
a=newa(rd,trd/colrad,colrad)
e=newe(rd,trd/colrad,colrad)
print('newa /au',a/au)
print('newe',e)
'''
