import numpy as np
from numpy import linalg as nplg
import matplotlib.pyplot as plt
from KCollisions.NewOrbit import *
import timeit
from numba import jit, int32,float64
import cProfile
import re
import csv

savename='DistTempSave.npy'
#Units
Mstar = 1.2e30
au = 1.496e11
year = 3.154e+7
#Setup
x='a'
a1=2*au
e1=0.99
s1=0
a2=2.1*au
e2=0.993
s2=1
dfSize = 1  #1 meter

#Rbeam = 2.6e3 #from paper veras
#TotalVol=2.26e14
Rbeam=1e3
TotalVol=1e10


@jit
def RcolCalc(N):
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
    V1 = np.array([rd1, R * td1])
    v1 = nplg.norm(V1)
    rd2 = rdot(a2, e2, C - s2, Mstar)
    td2 = thetadot(a2, e2, C - s2, Mstar)
    V2 = np.array([rd2, R * td2])
    v2 = nplg.norm(V2)
    Vrel = V2 - V1
    vrel=nplg.norm(Vrel)

    # Creating the TEmporary Momentum Data
    TempMomData = np.zeros((5, N))  # 0-rd,1-R*td,2-v,3-alpha,4-eccentricity
    H = np.linspace(0, 1, N)

    for i in range(0, N):
        TempMomData[0:2, i] = V1 + H[i] * Vrel # rd,R*td
        TempMomData[2, i] = nplg.norm(TempMomData[0:2, i])  # v
        TempMomData[3,i]=np.arccos(np.vdot(TempMomData[0:2, i],Vrel)/(TempMomData[2, i]*vrel))#alphacalc
        TempMomData[4,i]=newe(TempMomData[0, i],TempMomData[1, i]/R,R,Mstar)

    # Rcol
    RcolMat = np.zeros((N, N))
    PConstant = (3 * ((G * Mstar) ** 4.) * pi) / (256 * Rbeam * (R ** 12.))
    for i in range(0,N):#range(0,N-1):
        for j in range(0,N):#range(i+1,N):
            RcolMat[i,j]=(((1-TempMomData[4,i])**1.5)*((1-TempMomData[4,j])**1.5))/(TempMomData[2,i]*abs(sin(TempMomData[3,i]))*(((TempMomData[1,i]/R)**3.)*(TempMomData[1,j]/R)**3.))

    RcolMat=RcolMat*PConstant
    return RcolMat

@jit
def InitialDistCalc(J):
    IniDist=np.zeros((J+1))
    for f in range(1, J + 1):
        IniDist[f] = (2 * f * dfSize) ** (-3.5)
        IniDist[0] += IniDist[f]*(f**3.)
    #print('TotalVol unscaled',Dist[0,0,0])
    IniDist[1:J+1] = ((TotalVol / IniDist[0])) * IniDist[1:J + 1]
    IniDist[0]=TotalVol
    return IniDist
@jit
def ncolMatCalc(N,J):
    MassFraction = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            MassFraction[i, j] = (i ** 3.) / ((i ** 3.) + (j ** 3.))

    # collision bin calc
    ncolmat = np.zeros((N, N, J + 1, J + 1))
    '''
    ForcedAlt
    for n1 in range(0, N - 1):  # fragment 1
        for n2 in range(n1 + 1, N):
            ncolmat[n1, n2, :, :] = n1 + np.around(MassFraction[:, :] * (n2 - n1))
            for i in range(1, J + 1):
                ncolmat[n1, n2, i, i] = n1 + np.around(
                    ((n2 - n1) / 2) + ((-1) ** (np.floor(n2 / 2) + np.floor(n1 / 2))) * 0.1)
    '''
    for n1 in range(0, N - 1):  # fragment 1
        for n2 in range(n1 + 1, N):
            ncolmat[n1, n2, :, :] = n1 + np.around(MassFraction[:, :] * (n2 - n1))
    for n1 in range(0, N - 1):
        for i in range(1, J + 1):
            for n2 in range(n1 + 1, N - (n1)):
                ncolmat[n1, n2, i, i] = n1 + np.around((MassFraction[i, i] * (n2 - n1)) - 0.001)
                ncolmat[N - (n2 + 1), N - (n1 + 1), i, i] = N - (ncolmat[n1, n2, i, i] + 1)


    ncolmat = ncolmat.astype(int)
    #print(ncolmat[:,:,1,1])
    return ncolmat
@jit
def Evolution(N,J,Tpower):
    TimeStep = 1* year  # *(10**3.)
    T = 2*(10 ** Tpower)  # number of Timesteps
    #N = 500
    # Loading MomData
    RcolMat = RcolCalc(N)
    RcolMat = TimeStep * (dfSize ** 2.) *RcolMat #turning from a rate to a number of fragments

    # Setting up distributions
    TempT = int(10 ** 4.)
    print('TempT', TempT)
    ShowT = int(T / TempT)
    print('ShowT', ShowT)
    TempDist = np.zeros((J + 1, N, TempT + 1))
    Dist = np.zeros((J + 1, N, ShowT + 1))

    # Initial Parent Setup
    Dist[:,0,0]=InitialDistCalc(J)
    Dist[:, N - 1, 0] = Dist[:, 0, 0]

    # Setting up the collision bin results
    MassFraction = np.zeros((J + 1, J + 1))
    CrossSecMat = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            CrossSecMat[i, j] = (i + j) ** 2.


    ncolmat = ncolMatCalc(N,J)
     # Time Evolution
    for tshow in range(0, ShowT):  # tshow loop
        TempDist[1:(J+1), :, 0] = Dist[1:J+1, :, tshow]
        print('tshow', tshow)

        for t in range(0, TempT):  # temp loop
            TempDist[1:J+1, :, t + 1] = TempDist[1:J+1, :, t]  # new=old with changes
            # All Collision loop
            for n1 in range(0, N-1):  # velocity bin 1
                for n2 in range(n1+1, N):  # Velocity Bin 2
                    for f1 in range(1, J + 1):  # Fragment 1
                        for f2 in range(1, J + 1):  # Fragment 2
                            #print('n1,n2',n1,n2)
                            colnumber = RcolMat[n1,n2] * CrossSecMat[f2,f1]* TempDist[f1, n1, t] * TempDist[f2, n2, t]
                            # Add
                            TempDist[f2, ncolmat[n1,n2,f2,f1], t + 1] += colnumber
                            TempDist[f1, ncolmat[n1,n2,f2,f1], t + 1] += colnumber
                            #Remove
                            TempDist[f2, n2, t + 1] -= colnumber  # can get minus number of fragments, needs to remain super large
                            TempDist[f1, n1, t + 1] -= colnumber


        #looping back around
        Dist[1:J+1, :, tshow+1 ] = TempDist[1:J+1, :, TempT]
        # Total Mass Sum
        for f in range(1, J + 1):
            Dist[0, :, tshow + 1] += Dist[f, :, tshow + 1] * (f ** 3.)


    return Dist

@jit
def ncolMatCalcnoround(N,J):
    MassFraction = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            MassFraction[i, j] = (i ** 3.) / ((i ** 3.) + (j ** 3.))

    # collision bin calc
    ncolmattemp = np.zeros((N, N, J + 1, J + 1))
    for n1 in range(0, N - 1):  # fragment 1
        for n2 in range(n1 + 1, N):
            ncolmattemp[n1, n2, :, :] = n1 + MassFraction[:, :] * (n2 - n1)
    ncolmat=np.zeros((N,N,J+1,J+1,2))
    ncolmat[:,:,:,:,0]=np.floor(ncolmattemp)
    #ncolmat[:, :, :, :, 0]=ncolmat[:,:,:,:,0].astype(int)
    ncolmat[:, :, :, :, 1]=ncolmattemp-ncolmat[:,:,:,:,0]




    return ncolmat
@jit
def Evolutionbinsplit(N,J,Tpower):
    TimeStep = 1* year  # *(10**3.)
    T = 2*(10 ** Tpower)  # number of Timesteps
    #N = 500
    # Loading MomData
    RcolMat = RcolCalc(N)
    RcolMat = TimeStep * (dfSize ** 2.) *RcolMat #turning from a rate to a number of fragments

    # Setting up distributions
    TempT = int(10 ** 4.)
    print('TempT', TempT)
    ShowT = int(T / TempT)
    print('ShowT', ShowT)
    TempDist = np.zeros((J + 1, N, TempT + 1))
    Dist = np.zeros((J + 1, N, ShowT + 1))

    # Initial Parent Setup
    Dist[:,0,0]=InitialDistCalc(J)
    Dist[:, N - 1, 0] = Dist[:, 0, 0]

    # Setting up the collision bin results
    MassFraction = np.zeros((J + 1, J + 1))
    CrossSecMat = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            CrossSecMat[i, j] = (i + j) ** 2.


    ncolmat = ncolMatCalcnoround(N,J)
     # Time Evolution
    for tshow in range(0, ShowT):  # tshow loop
        TempDist[1:(J+1), :, 0] = Dist[1:J+1, :, tshow]
        print('tshow', tshow)

        for t in range(0, TempT):  # temp loop
            TempDist[1:J+1, :, t + 1] = TempDist[1:J+1, :, t]  # new=old with changes
            # All Collision loop
            for n1 in range(0, N-1):  # velocity bin 1
                for n2 in range(n1+1, N):  # Velocity Bin 2
                    for f1 in range(1, J + 1):  # Fragment 1
                        for f2 in range(1, J + 1):  # Fragment 2
                            #print('n1,n2',n1,n2)
                            colnumber = RcolMat[n1,n2] * CrossSecMat[f2,f1]* TempDist[f1, n1, t] * TempDist[f2, n2, t]
                            # Add
                            TempDist[f2, int(ncolmat[n1,n2,f2,f1,0]), t + 1] += (1-ncolmat[n1,n2,f2,f1,1])*colnumber
                            TempDist[f1, int(ncolmat[n1,n2,f2,f1,0]), t + 1] += (1-ncolmat[n1,n2,f2,f1,1])*colnumber

                            TempDist[f2, int(ncolmat[n1,n2,f2,f1,0])+1, t + 1] += colnumber*ncolmat[n1,n2,f2,f1,1]
                            TempDist[f1, int(ncolmat[n1,n2,f2,f1,0])+1, t + 1] += colnumber*ncolmat[n1,n2,f2,f1,1]

                            #Remove
                            TempDist[f2, n2, t + 1] -= colnumber  # can get minus number of fragments, needs to remain super large
                            TempDist[f1, n1, t + 1] -= colnumber


        #looping back around
        Dist[1:J+1, :, tshow+1 ] = TempDist[1:J+1, :, TempT]
        # Total Mass Sum
        for f in range(1, J + 1):
            Dist[0, :, tshow + 1] += Dist[f, :, tshow + 1] * (f ** 3.)


    return Dist
def savedistbinsplit(N,J,Tpower):
    Dist=Evolutionbinsplit(N,J,Tpower)
    np.save(savename,Dist)
    return



def savedist(N,J,Tpower):
    Dist=Evolution(N,J,Tpower)
    np.save(savename,Dist)
    return

def Timing():
    t = timeit.Timer(lambda: Evolution(100,4, 5.))
    print(t.timeit(number=1))
    return


#timeit.timeit(stmt=Timing(), number=1)
#ncolMattest=ncolMatCalc(11,1)
#print(ncolMattest[:,:,1,1])
#ncolMattest=ncolMatCalcnoround(11,1)
#print(ncolMattest[:,:,1,1,0])
#print(ncolMattest[:,:,1,1,1])
#savedistbinsplit(101,1,6.)
#savedist(201,1,6.)
#EvolutionGraphs()
#ParentMassGraphs()
#TotalMassGraphs()
#RcolPlot(500,10)

#plt.show()

#VrelPlot(500)
#Timing()

#TotalMassGraphs(500)
#Evolution(100,1,5)



#cProfile.runctx('Timing()',globals=globals(),locals=locals(),filename='Test.profile')
#python -m pstats ComponentSpaceAnalysis\Test.profile



