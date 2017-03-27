import matplotlib.pyplot as plt
import numpy as np

from KCollisions.KNewOrbit import *


def MinimalOrbitCascade(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,N): #Note at the moment N<=6, else the array is too large
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    if x=='a':
        R=CollisionData[1,0]
        C=CollisionData[0,0]
    elif x=='b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    #calculating the number of orbits in each generation
    n = np.zeros((2, N+1), dtype=int)  # number in row, total number     CHANGED
    n[0, 0] = 2
    n[1, 0] = 2

    for i in range(1, N+1):
        n[0, i] = 2**(i-1)
        n[1, i] = n[1, i - 1] + n[0, i]
    print('Total Number of Orbits'), print(n[1,N])
    # datasetup
    if n[0,N]>2:
        LargestGenSize =int(n[0, N])
    else:
        LargestGenSize = 2
    CascadeTable = np.zeros((N + 1, LargestGenSize + 1, 7)) # generations vertically, different orbits horizontally, with each orbits data stored in the third dimension
    CascadeTable[:, 0, 0] = n[0, :]  # the first column stores how many orbits are in each generation
    CascadeTable[:, :, 4] = 99  # the default eccentricity before being over written for ease of finding the minimal
    CascadeTable[0, 1, :] = [1, 0, 0, a1, e1, s1,m1]  # The First parent orbits data (generation 0): Norbit,parent1,parent2,a,e,s,m
    CascadeTable[0, 2, :] = [2, 0, 0, a2, e2, s2, m2]  # The second parent orbits
    [a3, e3, s3, m3] = NewOrbit(a1, e1, s1, a2, e2, s2, Mstar, m1, m2,K, C, R)  # generating the initial created orbit
    CascadeTable[1, 1, :] = [3, 1, 2, a3, e3, s3, m3]  # Storing the first generation orbit
    Norbit = 4  # orbit label number
    # Cascade Loop
    # the first generation to be calculated in the loop is the second
    for i in range(2, N + 1):  # generation loop, starting at 1
        CascadeTable[i, 0, 0] =int(CascadeTable[i, 0, 0] )
        j = 1  # write position on line
        for k in range(1, 3):  # Only the two parent orbits!

            # the second parent loop

                for y in range(1, int(n[0, i-1])+1 ):  # loops horizontally through the row for the second parent

                    [a5, e5, s5, m5] = CascadeTable[i-1, y, 3:7]  # the orbit data of the second parent

                    [a3, e3, s3, m3] = NewOrbit(CascadeTable[0,k,3],CascadeTable[0,k,4],CascadeTable[0,k,5], a5, e5, s5, Mstar, CascadeTable[0,k,6], m5,K, C,R)  # calculating the generated orbit
                    CascadeTable[i, j, :] = [Norbit, k, CascadeTable[i-1, y, 0], a3, e3, s3, m3]
                    j += 1  # moving the output along horizontally to write the next in the generation
                    Norbit += 1  # Labeling the orbit





    return CascadeTable  # Output is the entire data table

def MinimalCircularisation(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,N): #function finds the most circular orbits of each generation
    CascadeTable=MinimalOrbitCascade(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,N) #creating and loading the data table
    MinGen=np.zeros((N+1,8))
    for i in range(0,N+1): #for every generation, find the orbit with the min eccentricity and store the data, including the horizontal index for future reference
        MinIndex=CascadeTable[i,:,4].argmin()
        MinGen[i,:]=np.append(MinIndex, CascadeTable[i, MinIndex, :])

    return MinGen



def MinimalCircularisationPolarGraph(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,N):
    MinGen=MinimalCircularisation(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,N)
    F = np.linspace(0, 2 * pi, 100)
    plt.figure()
    #Plotting Parent Orbits
   # R = npr(a1, e1, F - s1)
    #plt.polar(F, R,label="Orbit1")
    #R = npr(a2, e2, F - s2)
    #plt.polar(F, R, label="Orbit2")

    for i in range(1,N+1):
        [a3,e3,s3]=MinGen[i,4:7]
        R = npr(a3, e3, F - s3)
        plt.polar(F, R, label='Generation %s'%(i))
    plt.legend()
    plt.title('Miminal Polar, with K=%s' %K)
    plt.show()
    return

def MinimalCircularisationEccentricityGraph(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,N):
    MinGen=MinimalCircularisation(x,a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,N)
    Generation = np.linspace(0, N, N+1)
    plt.figure()
    plt.plot(Generation,MinGen[:,5])
    plt.xlabel('Generation of Orbits')
    plt.ylabel('Eccentricity of the most circular orbit')
    plt.title('Minimal -Lowest Eccentricity of the Generation with K=%s'%K)
    plt.show()
    return


def MinimalEccentricitySGraph(x,a1,e1,a2,e2,Mstar,m1,m2,K,N,Ns):
    s1=0
    Generation = np.linspace(0, N, N + 1)
    plt.figure()
    steps=(2*pi)/Ns
    for s2 in range(1,Ns):

        MinGen = MinimalCircularisation(x,a1, e1, s1, a2, e2, np.around((steps*s2+0.001), decimals=2), Mstar, m1, m2,K, N)
        plt.plot(Generation, MinGen[:, 5],label='Lambda/pi=%s'%(np.around((steps*s2+0.001)/pi, decimals=2)))

    plt.xlabel('Generation of Orbits')
    plt.ylabel('e of the most circular orbit')
    plt.title('Lowest e of the Generation with K=%s, collision %s'%(K,x))
    plt.legend()
    plt.xlim([0,N])
    plt.show()
    return
