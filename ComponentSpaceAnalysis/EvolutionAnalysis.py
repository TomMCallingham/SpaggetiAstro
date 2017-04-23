from ComponentSpaceAnalysis.Evolution import  *
savename='DistTempSave.npy'
def ParentMassGraphs():
    Dist = np.load(savename)
    N = np.size(Dist[0, :, 0])
    ShowT = np.size(Dist[0, 0, :]) - 1
    TimeStep=year
    ShowTspace = np.linspace(0, TimeStep*1e6, ShowT + 1)
    Rcollog = 20
    testfunc = (2 * TotalVol) / (1 + np.exp(2 * TotalVol * ShowTspace * (10 ** -Rcollog)))
    plt.figure()
    plt.semilogy(ShowTspace, Dist[0, 0, :], label='First Parent Mass')
    plt.semilogy(ShowTspace, Dist[0, N - 1, :], label='Second Parent Mass')
    #plt.semilogy(ShowTspace,(testfunc), label='mass guess')
    plt.legend()
    '''
    plt.figure()
    plt.plot(ShowTspace, Dist[0, 0, :], label='First Parent Mass')
    plt.plot(ShowTspace, Dist[0, N - 1, :], label='Second Parent Mass')
    plt.plot(ShowTspace, (testfunc), label='mass guess')
    plt.legend()
    '''
    plt.show()
    return

def TotalMassGraphs():
    Dist = np.load(savename)
    N=np.size(Dist[0,:,0])
    ShowT = np.size(Dist[0, 0, :])-1
    TotalMass=np.zeros(ShowT+1)
    for t in range(0,ShowT+1):
        TotalMass[t]= np.sum(Dist[0,:,t])
    TotalMassChange=TotalMass-TotalMass[0]
    TimeStep=year
    ShowTspace = np.arange(0,  ShowT+1 )
    plt.figure()
    plt.plot(ShowTspace, TotalMassChange, label='Total Mass')

    plt.legend()
    plt.show()
    return

def EvolutionGraphs():
    Dist = np.load(savename)
    N = np.size(Dist[0, :, 0])
    ShowT = np.size(Dist[0, 0, :])
    H = np.linspace(0, 1, N)
    #plt.figure()
    Graphnum=10
    step=int(round(ShowT/Graphnum))
    for t in range(0, ShowT,step):
        plt.figure()
        plt.plot(H, Dist[0, :, t], label='t=%se5years' % t)
        plt.xlabel('H')
        plt.ylabel('TotalVolume')
        plt.legend()
        plt.ylim([0, 1.5 * TotalVol])

    plt.show()
    return



def RcolPlot(N,GraphNumber):
    RcolMat=RcolCalc(N)
    H = np.linspace(0, 1, N)
    plt.figure()
    step=int(N/GraphNumber)
    for i in range(0,GraphNumber-1):
        #plt.plot(H[(i*step)+1:N],np.log10(RcolMat[(i*step),(i*step)+1:N]),label='Rcol of H=%s'%(i*step/N))
        plt.plot(H, np.log10(RcolMat[(i * step), :]),label='Rcol of H=%s' % (i * step / N))
    plt.xlabel('H')
    plt.ylabel('Rcol/(ri^2+rj^2')
    plt.legend()
    plt.show()
    return

def EccentricityGraph(l):
    N=500
    # CollisionPoints
    s2=l
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    # setting up the parent data
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

    # Creating the TEmporary Momentum Data
    TempMomData = np.zeros((3, N))  # 0-rd,1-R*td,2-eccentricity
    H = np.linspace(0, 1, N)

    for i in range(0, N):
        TempMomData[0:2, i] = V1 + H[i] * Vrel # rd,R*td
        TempMomData[2,i]=newe(TempMomData[0,i],TempMomData[1,i]/R,R,Mstar)
    plt.plot(H,TempMomData[2,:],label='l=%s'%l)
    plt.legend()
    #plt.show()

    return
def eccentricitygraph():
    Dist = np.load(savename)
    N = np.size(Dist[0, :, 0])
    ShowT = np.size(Dist[0, 0, :])
    H = np.linspace(0, 1, N)

    return

def TotalMomentumPlot():
    Dist = np.load(savename)
    N = np.size(Dist[0, :, 0])
    ShowT = np.size(Dist[0, 0, :])
    H = np.linspace(0, 1, N)
    MomCentre=np.zeros(ShowT)
    for t in range(0,ShowT):
        for n in range(0,N):
            MomCentre[t] += Dist[0,n,t]*n
    MomCentre=MomCentre/(2*TotalVol)
    ShowTspace = np.arange(0, ShowT)
    plt.figure()
    plt.plot(MomCentre,ShowTspace,label='MemCentre')
    plt.show()
    return

def DifDistSave(Dist1name,Dist2name):
    Dist1=np.load(Dist1name)
    Dist2 = np.load(Dist2name)
    return

#eccentricitygraphs()

#EvolutionGraphs()
#ParentMassGraphs()
#TotalMassGraphs()
RcolPlot(500,10)
#TotalMomentumPlot()