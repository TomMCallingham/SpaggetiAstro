def BestCircularisation(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N):
    CascadeTable=OrbitCascade(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N)

    MinGen=np.zeros((2,N+1))

    for i in range(0,N+1):
        MinGen[0,i]=np.min(CascadeTable[i,:,4])
        MinGen[1,i]=CascadeTable[i,:,4].argmin()

    BestOrbitGen=MinGen[0,:].argmin()
    BestOrbitIndex=int(MinGen[1,BestOrbitGen])
    print('BestOrbitGen'),print(BestOrbitGen)
    BestOrbit=np.append(BestOrbitGen,CascadeTable[BestOrbitGen,BestOrbitIndex,:])
    print(BestOrbit)
    print('Mingen='),print(MinGen)
    return BestOrbit

def BestCircularistationGraph(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N):
    [Gen,Norbit,parent1,parent2,a3, e3, s3, m3]=BestCircularisation(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N)
    F = np.linspace(0, 2 * pi, 100);

    plt.figure()

#Plotting1
    R = npr(a1, e1, F - s1)
    plt.polar(F, R,label="Orbit1")
    R = npr(a2, e2, F - s2)
    plt.polar(F, R, label="Orbit2")
    R = npr(a3, e3, F - s3)
    plt.polar(F, R, label="BestOrbit")
    plt.legend()
    plt.show(1)
    return

def OldlessdataCircularisation(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N): #function finds the most circular orbits of each generation
    CascadeTable=OrbitCascade(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,N) #creating and loading the data table
    MinGen=np.zeros((2,N+1))
    for i in range(0,N+1): #for every generation
        MinGen[0,i]=np.min(CascadeTable[i,:,4]) #the min eccentricity of the ith generation
        MinGen[1,i]=CascadeTable[i,:,4].argmin() #the location of that orbit

    return MinGen
