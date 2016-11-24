from NewOrbit import *
import matplotlib.pyplot as plt
import numpy as np

def OrbitCascade(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N): #Note at the moment N<=6, else the array is too large
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    #for now, just selecting closest collision point
    if CollisionData[1,0]>CollisionData[1,1]:
        R=CollisionData[1,1]
        C=CollisionData[0,1]
    else:
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]

    #calculating the number of orbits in each generation
    n = np.zeros((2, N+1))#, dtype=int)  # number in row, total number     CHANGED
    n[0, 0] = 2
    n[0, 1] = 1
    n[0, 2] = 2
    n[1, 0] = 2
    n[1, 1] = 3
    n[1, 2] = 5
    i = 3
    for i in range(3, N+1):
        n[0, i] = n[0, i - 1] * n[1, i - 2]
        n[1, i] = n[1, i - 1] + n[0, i]
        i += 1
    print('Total Number of Orbits'), print(n[1,N])
    #datasetup
    CascadeTable=np.zeros((N+1,int(n[0,N])+1,7)) #generations vertically, different orbits horizontally, with each orbits data stored in the third dimension
    CascadeTable[:,0,0]=n[0,:] #the first column stores how many orbits are in each generation
    CascadeTable[:, :, 4]=99 #the default eccentricity before being over written for ease of finding the minimal
    CascadeTable[0,1,:]=[1,0,0,a1,e1,s1,m1] #The First parent orbits data (generation 0): Norbit,parent1,parent2,a,e,s,m
    CascadeTable[0,2,:]=[2,0,0,a2,e2,s2,m2]#The second parent orbits
    [a3, e3, s3,m3]=NewOrbit(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,C,R) # generating the initial created orbit
    CascadeTable[1,1,:]=[3,1,2,a3, e3, s3,m3] #Storing the first generation orbit
    Norbit=4  # orbit label number
    #Cascade Loop
    i=2 # the first generation to be calculated in the loop is the second
    for i in range(2, N+1): #generation loop
        k=1 #first parent position in line above
        j=1 #write position on line
        for k in range(1,int(n[0,i-1])+1):  #first parent loop, from the line above
            [a4,e4,s4,m4]=CascadeTable[i-1,k,3:7] #fixing the first parent from line above
            #the second parent loop
            x=0 # x represents the generation of the second parent
            for x in range(0,i-1): #stops before i! looping down vertically for second parent
                y=1 #y reperesents the horizontal placement of the second parent
                for y in range(1,int(n[0,x])+1): #loops through that row for the second parent
                    [a5,e5,s5,m5]=CascadeTable[x,y,3:7] #the orbit data of the second parent
                    [a3, e3, s3,m3] = NewOrbit(a4, e4, s4, a5, e5, s5, Mstar, m4, m5, C, R) #calculating the generated orbit
                    CascadeTable[i, j, :] = [Norbit, CascadeTable[i-1,k,0], CascadeTable[x,y,0], a3, e3, s3,m3]
                    j += 1 #moving the output along horizontally to write the next in the generation
                    Norbit += 1 #Labeling the orbit
                    y +=1 #selecting the next second parent horizontally
                x +=1 #selecting the next generation of parents
            k += 1 #selecting the next first parent
    return CascadeTable #Output is the entire data table

def Circularisation(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N): #function finds the most circular orbits of each generation
    CascadeTable=OrbitCascade(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N) #creating and loading the data table
    MinGen=np.zeros((N+1,8))
    for i in range(0,N+1): #for every generation, find the orbit with the min eccentricity and store the data, including the horizontal index for future reference
        MinIndex=CascadeTable[i,:,4].argmin()
        MinGen[i,:]=np.append(MinIndex, CascadeTable[i, MinIndex, :])

    return MinGen



def CircularisationPolarGraph(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N):
    MinGen=Circularisation(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N)
    F = np.linspace(0, 2 * pi, 100)
    plt.figure()
    #Plotting Parent Orbits
    R = npr(a1, e1, F - s1)
    plt.polar(F, R,label="Orbit1")
    R = npr(a2, e2, F - s2)
    plt.polar(F, R, label="Orbit2")
    i=1
    for i in range(1,N+1):
        [a3,e3,s3]=MinGen[i,4:7]
        R = npr(a3, e3, F - s3)
        plt.polar(F, R, label='Generation %s'%(i))
    plt.legend()
    plt.show()
    return

def CircularisationEccentricityGraph(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N):
    MinGen=Circularisation(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N)
    Generation = np.linspace(0, N, N+1)
    plt.figure()
    plt.plot(Generation,MinGen[:,5])
    plt.xlabel('Generation of Orbits')
    plt.ylabel('Eccentricity of the most circular orbit')
    plt.title('Lowest Eccentricity of the Generation')
    plt.show()
    return




