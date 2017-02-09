import numpy as np
from numpy import linalg as nplg
import matplotlib.pyplot as plt
from KCollisions.NewOrbit import *
import timeit
from numba import jit, int32,float64
import cProfile
import re


#Units
Mstar = 1.2e30
au = 1.496e11
year = 3.154e+7
#Setting up collisions!


x='a'
a1=2*au
e1=0.99
s1=0
a2=2.1*au
e2=0.993
s2=2
dfSize = 1  # 1/(1.496e11) #1 meter
Rbeam = 1000 * dfSize

def MomDataCalc( N):


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
    plt.show()
    '''

    return MomData


@jit
def Evolution(N):

    TimeStep = 0.1* year  # *(10**3.)
    T = 10 ** 6  # number of Timesteps
    J = 4  # fragment radius bins
    #N = 500
    # Loading MomData
    MomData = MomDataCalc(N)
    MomData[5:7, :] = TimeStep * MomData[5:7, :]

    # Setting up distributions
    TempT = int(10 ** 4.)
    print('TempT', TempT)
    ShowT = int(T / TempT)
    print('ShowT', ShowT)
    TempDist = np.zeros((J + 1, N, TempT + 1))
    Dist = np.zeros((J + 1, N, ShowT + 1))

    # Parent Setup
    # Starting fragment distribution
    TotalMass=10**6.
    for f in range(1, J + 1):
        Dist[f, 0, 0] = (2 * f * dfSize) ** (-3.5)
        Dist[0, 0, 0] += Dist[f, 0, 0]*(f**3.)
    #print('TotalMass unscaled',Dist[0,0,0])
    Dist[1:J+1, 0, 0] = ((TotalMass/Dist[0,0,0])) * Dist[1:J+1, 0, 0]
    Dist[0, 0, 0]=10**6.
    Dist[:, N - 1, 0] = Dist[:, 0, 0]
    #
    print('largest number particles:',np.max(Dist[1:N,0,0]))
    '''print('Initial Dist',Dist[:,0,0])'''

    nparent = [0,N-1]


    # Settign up the collision bin results
    MassFraction = np.zeros((J + 1, J + 1))
    CrossSecMat = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            MassFraction[i, j] = (i ** 3.) / ((i ** 3.) + (j ** 3.))
            CrossSecMat[i, j] = (i + j) ** 2.

    #parent collisions
    ncolmatparent = np.zeros((J + 1, J + 1))
    ncolmatparent = np.round((N - 1) * MassFraction)
    ncolmatparent = ncolmatparent.astype(int)
    '''print('CrossSecMat',CrossSecMat)'''
    #child collisions
    ncolmat=np.zeros((2,N,J+1,J+1))
    for p in (0, 1):  # parents 1 or 2
        for n in range(1, N - 1):  # Velocity Bin Loop
            ncolmat[p,n,:,:] = n + np.round(MassFraction * (nparent[p] - n))
    ncolmat = ncolmat.astype(int)
    #setting up rcoll
    #Rcolmat = np.zeros((J + 1, J + 1))
    Rcolmat = MomData[5, N - 1] * 0.5*(np.max(Dist[1:N,0,0])**2.)
    print('parent MomData',MomData[5,N-1])
    print('first collisoin',Rcolmat)

    # Time Evolution
    for tshow in range(0, ShowT):  # tshow loop

        TempDist[1:(J+1), :, 0] = Dist[1:J+1, :, tshow]
        print('tshow', tshow)

        for t in range(0, TempT):  # temp loop
            TempDist[1:J+1, :, t + 1] = TempDist[1:J+1, :, t]  # new=old with changes


            #Parent-Parent Collision Loop
            for fp in range(1, J + 1):
                for f in range(1, J + 1):
                    Rcolmat = MomData[5, N - 1] * CrossSecMat[f,fp] * TempDist[f, 0, t] * TempDist[fp, N - 1, t]
                    # Add
                    TempDist[f, ncolmatparent[fp, f], t + 1] += Rcolmat
                    TempDist[fp, ncolmatparent[fp, f], t + 1] += Rcolmat
                    # Remove
                    TempDist[f, 0, t + 1] -= Rcolmat
                    TempDist[fp, N - 1, t + 1] -= Rcolmat

            '''
            # Parent-Child Collision loop
            for p in (0, 1):  # parents 1 or 2
                npar=nparent[p]
                for n in range(1, N - 1):  # Velocity Bin Loop
                    for fp in range(1, J + 1):  # Parent Fragment loop
                        for f in range(1, J + 1):  # ChildFragment Loop
                            Rcolmat = MomData[5 + p, n] * CrossSecMat[f,fp]* TempDist[f, n, t] * TempDist[fp, npar, t]
                            # Add
                            TempDist[f, ncolmat[p,n,fp,f], t + 1] += Rcolmat
                            TempDist[fp, ncolmat[p,n,fp,f], t + 1] += Rcolmat
                            #Remove
                            TempDist[f, n, t + 1] -= Rcolmat  # can get minus number of fragments, needs to remain super large
                            TempDist[fp, npar, t + 1] -= Rcolmat

            '''




        #looping back around
        Dist[1:J+1, :, tshow+1 ] = TempDist[1:J+1, :, TempT]
        # Total Mass Sum
        '''Dist[0,:,tshow+1]=0'''
        for f in range(1, J + 1):
            Dist[0, :, tshow + 1] += Dist[f, :, tshow + 1] * (f ** 3.)


    return Dist

def ParentMassGraphs(N):
    Dist = Evolution(N)
    ShowT = np.size(Dist[0, 0, :])
    TimeStep=year
    ShowTspace = np.linspace(0, TimeStep, ShowT )#+ 1)
    plt.figure()
    plt.plot(ShowTspace, Dist[0, 0, :], label='First Parent Mass')
    plt.plot(ShowTspace, Dist[0, N - 1, :], label='Second Parent Mass')
    plt.legend()
    plt.show()
    return
def TotalMassGraphs(N):
    Dist = Evolution(N)
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

def EvolutionGraphs(N):
    Dist = Evolution(N)
    ShowT = np.size(Dist[0, 0, :])
    H = np.linspace(0, 1, N)
    #plt.figure()
    Graphnum=5
    step=int(round(ShowT/Graphnum))
    for t in range(0, ShowT,step):
        plt.figure()
        plt.plot(H, Dist[0, :, t], label='t=%s' % t)
        plt.legend()

    plt.show()
    return


def Timing():
    N=500
    Evolution(N)
    return




#EvolutionGraphs(500)
#ParentMassGraphs(500)
TotalMassGraphs(500)


#print(timeit.timeit(Timing, number=3))
#cProfile.runctx('Timing()',globals=globals(),locals=locals(),filename='Test.profile')
#python -m pstats ComponentSpaceAnalysis\Test.profile
