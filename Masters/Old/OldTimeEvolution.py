import numpy as np
from numpy import linalg as nplg
import matplotlib.pyplot as plt
from KCollisions.KNewOrbit import *
import timeit

Mstar = 1.2e30
au = 1.496e11
year = 3.154e+7


def MomDataCalc(x, a1, e1, s1, a2, e2, s2, N):
    dfSize = 1  # 1/(1.496e11) #1 meter
    Rbeam = 1000 * dfSize

    # CollisionPoints
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
    V1 = np.array([[rd1], [R * td1]])
    v1 = nplg.norm(V1)
    rd2 = rdot(a2, e2, C - s2, Mstar)
    td2 = thetadot(a2, e2, C - s2, Mstar)
    V2 = np.array([[rd2], [R * td2]])
    v2 = nplg.norm(V2)
    Vrel = V2 - V1

    # Creating the TEmporary Momentum Data
    TempMomData = np.zeros((8, N))  # 0-rd,1-R*td,2-v,3-vrel1,4-vrel2,5-angle1,6-angle2,7-e,
    H = np.linspace(0, 1, N)
    TempMomData[0:2, :] = V1 + H * Vrel  # rd,td

    # Finding The Collision Angles, skipping parents
    for i in range(0, N):
        TempMomData[2, i] = nplg.norm(TempMomData[0:2, i])  # v
        TempMomData[3, i] = nplg.norm(TempMomData[0:2, i] - V1)  # vrel1
        TempMomData[4, i] = nplg.norm(TempMomData[0:2, i] - V2)  # vrel2
        TempMomData[5, i] = np.arccos(np.vdot(TempMomData[0:2, i], V1) / (v1 * TempMomData[2, i]))  # angle 1
        TempMomData[6, i] = np.arccos(np.vdot(TempMomData[0:2, i], V2) / (v2 * TempMomData[2, i]))  # angle2
    # e calc
    TempMomData[7, :] = newe(TempMomData[0, :], TempMomData[1, :] / R, R, Mstar)  # eccentricity
    # Rcol
    RcolParents = np.zeros((2, N))
    MomentaCol = np.zeros((2, N))  # orbit dependand part, 1 for each parent
    PConstant = (3 * ((G * Mstar) ** 4.) * pi) / (256 * Rbeam * (R ** 12.))
    ParentConstant = np.zeros((2))
    ParentConstant[0] = 1 / (v1 * ((td1 ** 3.)))  # parent1 constant
    ParentConstant[1] = 1 / (v2 * ((td2 ** 3.)))  # parent1 constant
    # Note parents cannot collide with themselvels
    MomentaCol[0, 1:N] = (TempMomData[3, 1:N] / (
    TempMomData[2, 1:N] * ((TempMomData[1, 1:N] / R) ** 3.) * abs(np.sin(TempMomData[5, 1:N]))))  # parent1
    MomentaCol[1, 0:N - 1] = (TempMomData[4, 0:N - 1] / (
    TempMomData[2, 0:N - 1] * ((TempMomData[1, 0:N - 1] / R) ** 3.) * abs(np.sin(TempMomData[6, 0:N - 1]))))  # paernt2
    RcolParents[0, :] = PConstant * ParentConstant[0] * MomentaCol[0, :] * (dfSize ** 2.)
    RcolParents[1, :] = PConstant * ParentConstant[1] * MomentaCol[1, :] * (dfSize ** 2.)

    # Creating the Output Data
    MomData = np.zeros((7, N))  # 0-rd,1-td,2-e,3-vrel1,4-vrel2,5-Rcol1,6Rcol2
    MomData[0:2, :] = TempMomData[0:2, :]  # momentums
    MomData[2, :] = TempMomData[7, :]  # eccentricity
    MomData[3:5, :] = TempMomData[3:5, :]  # Vrels
    MomData[5:7, :] = RcolParents  # Col prob without the radius squared in dfsize units
    '''
    plt.figure()
    plt.plot(H[1:N],np.log10(MomData[5,1:N]),label='Rcol parent1')
    plt.plot(H[0:N-1], np.log10(MomData[6, 0:N-1]), label='Rcolparent 2')
    plt.legend()
    '''
    return MomData


def Evolution(x, a1, e1, s1, a2, e2, s2):
    TimeStep = 10 * year  # *(10**3.)
    T = 10 ** 4  # number of Timesteps
    J = 4  # fragment radius bins
    N = 100
    dfSize = 1
    # Loading MomData
    MomData = MomDataCalc(x, a1, e1, s1, a2, e2, s2, N)
    MomData[5:7, :] = TimeStep * MomData[5:7, :]

    # Setting up distributions
    TempT = int(10 ** 3.)
    print('TempT', TempT)
    ShowT = int(T / TempT)
    print('ShowT', ShowT)
    TempDist = np.zeros((J + 1, N, TempT + 1))
    Dist = np.zeros((J + 1, N, ShowT + 1))

    # Parent Setup
    # Starting fragment distribution
    for i in range(1, J + 1):
        Dist[i, 0, 0] = (2 * i * dfSize) ** (-3.5)
        Dist[0, 0, 0] += Dist[i, 0, 0]

    Dist[:, 0, 0] = (10 ** 6.) * Dist[:, 0, 0]
    Dist[:, N - 1, 0] = Dist[:, 0, 0]
    nparent = np.array([0, N - 1])  # N - 1])
    # Radius Calcs
    MassFraction = np.zeros((J + 1, J + 1))
    CrossSec = np.zeros((J + 1, J + 1))

    for i in range(1, J + 1):
        for j in range(1, J + 1):
            MassFraction[i, j] = (i ** 3.) / ((i ** 3.) + (j ** 3.))
            CrossSec[i, j] = (i + j) ** 2.


    ncolmat = np.zeros((J + 1, J + 1))  # ,dtype=int)
    Rcolmat = np.zeros((J + 1, J + 1))
    ncolmatparent = np.round((N - 1) * MassFraction)
    ncolmatparent = ncolmat.astype(int)
    ncolmat=np.zeros((2,N,J+1,J+1))
    for p in (0, 1):  # parents 1 or 2
        for n in range(1, N - 1):  # Velocity Bin Loop
            ncolmat[p,n,:,:] = n + np.round(MassFraction * (nparent[p] - n))
    ncolmat = ncolmat.astype(int)



    # Time Evolution
    for tshow in range(0, ShowT):  # tshow loop
        TempDist[:, :, 0] = Dist[:, :, tshow]
        print('tshow', tshow)

        for t in range(0, TempT):  # temp loop
            TempDist[:, :, t + 1] = TempDist[:, :, t]  # new=old with changes

            # Parent-Parent Collision Loop

            Rcolmat=MomData[5,N-1]*CrossSec*np.outer(TempDist[:,0,t],TempDist[:,N-1,t])
            #print('ncolmatunrounded',ncolmatunrounded)
            #print('ncolmatunrounded[1,1]',ncolmatunrounded[1,1])
            for fp in range(1, J + 1):
                for f in range(1, J + 1):
                    # find output bin
                    '''ncol = int(round(MassFraction[fp, f] * (N - 1)))
                    # ColRate
                    Rcol = MomData[5, N - 1] * CrossSec[f, fp] * TempDist[f, 0, t] * TempDist[fp, N - 1, t]'''

                    # changes

                    TempDist[f, ncolmatparent[fp,f], t + 1] +=  Rcolmat[f,fp]
                    TempDist[fp, ncolmatparent[fp,f], t + 1] +=  Rcolmat[f,fp]

                    TempDist[f, 0, t + 1] -= Rcolmat[f,fp]  # can get minus number of fragments, needs to remain super large
                    TempDist[fp, N - 1, t + 1] -= Rcolmat[f,fp]

            # Parent-Child Collision loop
            for p in (0, 1):  # parents 1 or 2
                for n in range(1, N - 1):  # Velocity Bin Loop
                    '''ncolmatunrounded=n +np.round(MassFraction * (nparent[p] - n))
                    ncolmatunrounded = ncolmatunrounded.astype(int)'''
                    Rcolmat = MomData[5, N - 1] * CrossSec * np.outer(TempDist[:, n, t], TempDist[:, nparent[p], t])
                    for fp in range(1, J + 1):  # Parent Fragment loop
                        for f in range(1, J + 1):  # ChildFragment Loop
                            # find output bin
                            '''ncol = int(n + round(MassFraction[fp, f] * (nparent[p] - n)))
                            # ColRate
                            Rcol = MomData[5 + p, n] * CrossSec[f, fp] * TempDist[f, n, t] * TempDist[fp, nparent[p], t]'''

                            # changes
                            TempDist[f, ncolmat[p,n,fp,f], t + 1] += Rcolmat[f,fp]
                            TempDist[fp, ncolmat[p,n,fp,f], t + 1] += Rcolmat[f,fp]

                            TempDist[
                                f, n, t + 1] -= Rcolmat[f,fp]  # can get minus number of fragments, needs to remain super large
                            TempDist[fp, nparent[p], t + 1] -= Rcolmat[f,fp]

            # Total Mass Sum
            for f in range(1, J + 1):
                TempDist[0, :, t + 1] = TempDist[0, :, t + 1] + TempDist[f, :, t + 1] * (f ** 3.)

        Dist[:, :, tshow + 1] = TempDist[:, :, TempT]

    return Dist


def EvolutionGraphs(x, a1, e1, s1, a2, e2, s2):
    Dist = Evolution(x, a1, e1, s1, a2, e2, s2)
    N = np.size(Dist[0, :, 0])
    ShowT = np.size(Dist[0, 0, :])
    H = np.linspace(0, 1, N)

    for t in range(0, ShowT):
        plt.figure()
        plt.plot(H, Dist[0, :, t], label='t=%s' % t)
        plt.legend()

    plt.show()
    return


def Timing():
    Evolution('a', 2 * au, 0.99, 0, 2.1 * au, 0.993, 2)
    return


print(timeit.timeit(Timing, number=1))
