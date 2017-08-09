'''
           #Vectorised
           TempDist[fp, ncolmatparent[fp, :], t + 1] += Rcolmat[:, fp]
           TempDist[fp, N - 1, t + 1] -= np.sum(Rcolmat[:, fp])

           TempDist[:, ncolmatparent[fp, :], t + 1] += Rcolmat[:, fp]
           TempDist[:, 0, t + 1] -= Rcolmat[:, fp]
           '''

# Parent-Parent Collision Loop
'''Rcolmat=MomData[5,N-1]*CrossSec*np.outer(TempDist[:,0,t],TempDist[:,N-1,t])
for fp in range(1, J + 1):

    for f in range(1, J + 1):
        #Add
        TempDist[f, ncolmatparent[fp,f], t + 1] +=  Rcolmat[f,fp]
        TempDist[fp, ncolmatparent[fp,f], t + 1] +=  Rcolmat[f,fp]
        #Remove
        TempDist[f, 0, t + 1] -= Rcolmat[f,fp]  # can get minus number of fragments, needs to remain super large
        TempDist[fp, N - 1, t + 1] -= Rcolmat[f,fp]
     # Parent-Child Collision loop
for p in (0, 1):  # parents 1 or 2
    for n in range(1, N - 1):  # Velocity Bin Loop
        Rcolmat = MomData[5+p, n] * CrossSec * np.outer(TempDist[:, n, t], TempDist[:, nparent[p], t])

        for fp in range(1, J + 1):  # Parent Fragment loop
            for f in range(1, J + 1):  # ChildFragment Loop

                # Add
                TempDist[f, ncolmatunrounded[p,n,fp,f], t + 1] += Rcolmat[f,fp]
                TempDist[fp, ncolmatunrounded[p,n,fp,f], t + 1] += Rcolmat[f,fp]
                #Remove
                TempDist[f, n, t + 1] -= Rcolmat[f,fp]  # can get minus number of fragments, needs to remain super large
                TempDist[fp, nparent[p], t + 1] -= Rcolmat[f,fp]
'''


#quicker even to calculate the cross sec fresh everytime!
'''Rcolmat = MomData[5, N - 1] * CrossSec[f, fp] * TempDist[f, 0, t] * TempDist[fp, N - 1, t]
CrossSec = np.zeros((J + 1, J + 1))
for i in range(1, J + 1):
    for j in range(1, J + 1):
        CrossSec[i, j] = (i + j) ** 2.'''

for fp in range(1, J + 1):
    for f in range(1, J + 1):
        speedfunc(MomData[5, N - 1], CrossSec[f, fp], TempDist[f, 0, t], TempDist[fp, N - 1, t])
        Rcolmat = MomData[5, N - 1] * CrossSec[f, fp] * TempDist[f, 0, t] * TempDist[fp, N - 1, t]
        '''Rcolmat = MomData[5, N - 1] * ((f+fp)**2.) * TempDist[f, 0, t] * TempDist[fp, N - 1, t]'''
        # Add
        TempDist[f, ncolmatparent[fp, f], t + 1] += Rcolmat
        TempDist[fp, ncolmatparent[fp, f], t + 1] += Rcolmat
        # Remove
        TempDist[f, 0, t + 1] -= Rcolmat
        TempDist[fp, N - 1, t + 1] -= Rcolmat
        # Parent-Child Collision loop



        '''
        def VrelPlot(N):
            MomData = MomDataCalc(N)
            H = np.linspace(0, 1, N)
            plt.figure()
            plt.plot(H[1:N], MomData[3, 1:N], label='Vrel parent1')
            plt.plot(H[0:N - 1], MomData[4, 0:N - 1], label='Vrel parent 2')
            plt.legend()
            #plt.show()
            return

              #Plot ncol test
            H = np.linspace(0, 1, N)
            plt.figure()
            GraphNumber=10
            step = int(N / GraphNumber)
            for i in range(0, GraphNumber - 1):
                plt.plot(H,ncolmatunrounded[i*step,:,1,1]/N ,label='ncol of H=%s' % (i * step / N))

            plt.legend()
            plt.show()
        '''

#Round Bins

@jit
def roundncolMatCalc(N,J):
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

def RoundEvolution(N,J,Tpower):
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