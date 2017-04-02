import numpy as np
from numpy import linalg as nplg
import matplotlib.pyplot as plt
from KCollisions.KNewOrbit import *
import timeit
from numba import jit, int32,float64
from Fragmentation.Fragmenting import *
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
Dmax=3 #poweres
Dmin=0
J=8
FragSizeDist = np.linspace(Dmin,Dmax,J)
FragSizeDist=10**FragSizeDist

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
    print('R in au',R)
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
    for i in range(0,N-1):
        for j in range(i+1,N):
            RcolMat[i,j]=(((1-TempMomData[4,i])**1.5)*((1-TempMomData[4,j])**1.5))/(TempMomData[2,i]*abs(sin(TempMomData[3,i]))*(((TempMomData[1,i]/R)**3.)*(TempMomData[1,j]/R)**3.))

    RcolMat=RcolMat*PConstant
    return RcolMat



@jit
def InitialDistCalc(J):
    IniDist=np.zeros((J+1))
    for f in range(1, J + 1):
        IniDist[f] = (FragSizeDist[f-1]) ** (-3.5)
        IniDist[0] += IniDist[f]*((FragSizeDist[f-1]/2)**3.)
    #print('TotalVol unscaled',Dist[0,0,0])

        IniDist = ((TotalVol / IniDist[0])) * IniDist[1:J + 1]
    return (IniDist)


@jit
def DistCalc(Botind,Topind):
    Dist = np.zeros((J + 1))
    for f in range(Botind, Topind+1):
        Dist[f] = (FragSizeDist[f - 1]) ** (-3.5)
        Dist[0] += Dist[f] * ((FragSizeDist[f - 1] / 2) ** 3.)
    # print('TotalVol unscaled',Dist[0,0,0])
    Dist[1:J + 1] = ((TotalVol / Dist[0])) * Dist[1:J + 1]
    Dist[0] = TotalVol
    return (Dist)

@jit
def ncolMatCalc(N, J):
    MassFraction = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            MassFraction[i, j] = (FragSizeDist[i-1] ** 3.) / ((FragSizeDist[i-1]  ** 3.) + (FragSizeDist[j-1]  ** 3.))

    # collision bin calc
    ncolmattemp = np.zeros((N, N, J + 1, J + 1))
    for n1 in range(0, N - 1):  # fragment 1
        for n2 in range(n1 + 1, N):
            ncolmattemp[n1, n2, :, :] = n1 + MassFraction[:, :] * (n2 - n1)
    ncolmat=np.zeros((N,N,J+1,J+1,2))
    ncolmat[:,:,:,:,0]=np.floor(ncolmattemp)  #lower bin
    #ncolmat[:, :, :, :, 0]=ncolmat[:,:,:,:,0].astype(int)
    ncolmat[:, :, :, :, 1]=ncolmattemp-ncolmat[:,:,:,:,0]  #fraction in lower bin
    return (ncolmat)

def colfragmatCalc(N,J):
    colfragmat=np.zeros((N,N,J+1,J+1,J+1))
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
    rd2 = rdot(a2, e2, C - s2, Mstar)
    td2 = thetadot(a2, e2, C - s2, Mstar)
    V2 = np.array([rd2, R * td2])
    Vrel = V2 - V1
    vrel = nplg.norm(Vrel)
    Qd=DispThres(FragSizeDist)
    for i in range(0,N-1):
        for j in range(i+1,N):
            Vcol=vrel*((j-i)/(N-1))
            for f1 in range(1, J ):  # Fragment 1
                for f2 in range(f1, J + 1):  # Fragment 2 larger





                    colfragmat[i,j,f1,f2,:]='test'
                    colfragmat[i, j, f2, f1, :]=colfragmat[i,j,f1,f2,:] #reflective



    return



@jit
def Evolution(N, J, TotalTime):

    #Load initialdist and Rcol
    RcolMat = RcolCalc(N)
    InitialDist = DistCalc(1, J) #Full Size bins
    InitialDist = ((TotalVol / InitialDist[0])) * InitialDist #Scale Initial Dist


    TimeStepMin=(1/(100*np.max(RcolMat)*2*(1**2.)*2*InitialDist[1]))/year
    TimeStepMinRound=int(10**floor(log10(TimeStepMin)))

    TimeStep = TimeStepMinRound * year
    #TotalTime=2 * (10 ** Tpower)
    T = int(TotalTime/TimeStepMinRound)  # number of Timesteps
    print('T',T)
    print('TimeStep',TimeStepMinRound)


    RcolMat = TimeStep *  RcolMat
    # Setting up distributions
    TempT = int(10 ** 3.)
    ShowT = int(T / TempT)
    print('TempT', TempT)
    print('ShowT', ShowT)
    print((TempT * TimeStep) / year)
    TempDist = np.zeros((J + 1, N, TempT + 1))
    Dist = np.zeros((J + 1, N, ShowT + 1))

    # Initial Parent Setup

    Dist[:,0,0]=InitialDist
    Dist[:, N - 1, 0] = Dist[:, 0, 0]

    # Setting up the collision bin results
    MassFraction = np.zeros((J + 1, J + 1))
    CrossSecMat = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            CrossSecMat[i, j] = (FragSizeDist[i-1] + FragSizeDist[j-1]) ** 2.


    ncolmat = ncolMatCalc(N, J)
    colfragmat=colfragmatCalc(N,J)
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
                            TempDist[:, int(ncolmat[n1,n2,f2,f1,0]), t + 1] += (1-ncolmat[n1,n2,f2,f1,1])*colnumber*colfragmat[n1,n2,f1,f2,:]

                            TempDist[:, int(ncolmat[n1, n2, f2, f1, 0])+1, t + 1] += colnumber*ncolmat[n1,n2,f2,f1,1]*colfragmat[n1,n2,f1,f2,:]


                            #Remove
                            TempDist[f2, n2, t + 1] -= colnumber  # can get minus number of fragments, needs to remain super large
                            TempDist[f1, n1, t + 1] -= colnumber


        #looping back around
        Dist[1:J+1, :, tshow+1 ] = TempDist[1:J+1, :, TempT]
        # Total Mass Sum
        for f in range(1, J + 1):
            Dist[0, :, tshow + 1] += Dist[f, :, tshow + 1] * ((dfSize*f) ** 3.)

    print((TempT * TimeStep)/year)
    return Dist




def savedist(N,J,Tpower):
    Dist=Evolution(N,J,Tpower)
    np.save(savename,Dist)
    return

def Timing():
    t = timeit.Timer(lambda: Evolution(100,4, 8.))
    print(t.timeit(number=1))
    return


RcolCalc(201)
#print(InitialDistCalc(1))
#timeit.timeit(stmt=savedist(201,10,2*(10**7.)), number=1)
#ncolMattest=ncolMatCalc(11,1)
#print(ncolMattest[:,:,1,1])
#ncolMattest=ncolMatCalcnoround(11,1)
#print(ncolMattest[:,:,1,1,0])
#print(ncolMattest[:,:,1,1,1])
#savedist(201,10,2*(10**7.))
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



