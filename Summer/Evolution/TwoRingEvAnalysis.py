from Summer.Evolution.TwoRingEvol import *

def MassChange(s1,s2,x,T,tstep,J):
    print('Evolving...')
    Dist=StaticEvolution(s1,s2,x,T,tstep,J)
    Tnumber = int(T / tstep)
    Time=np.arange(0,T+tstep,tstep)

    print('Plotting...')
    '''
    plt.figure()
    plt.plot(Time,Dist[0,0,:], label='P1 Vol')
    plt.plot(Time, Dist[1, 0, :], label='P2 Vol')
    plt.plot(Time, Dist[2, 0, :], label='C Vol')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Vol')
    '''
    plt.figure()
    plt.plot(Time,Dist[2,0,:]/(Dist[0,0,0]+Dist[1,0,0]),label='Vol Fraction')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Vol Fraction')
    '''
    #Sizes P1 Graph
    IniDist = InitialDistCalc(J)
    plt.figure()
    for f in range(1,J+1):
        plt.semilogx(Time,Dist[0,f,:]*(IniDist[1,f]**3.5),label='r=%s'%IniDist[1,f])
    plt.title('P1 Fragments over time')
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Vol')
    '''
    print('finshed')
    return
MassChange(0,3,'a',1e6,1,100)
plt.show()