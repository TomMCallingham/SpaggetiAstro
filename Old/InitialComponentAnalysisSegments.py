import matplotlib.pyplot as plt

from KCollisions.KNewOrbit import *

au= 1#149597871e3
G = 2.982e-27#6.67408e-11

def Segments(Indexes,Function): # finds the sections for which a graph is negative, given crossing points
    #print('Indexes',Indexes)

    NRanges=np.array([[0],[0]]) #creates the starting array, delete first after
    for i in range (0,np.size(Indexes)):  #note doesnt do the last, not that it matters
        if (Function[Indexes[i]])>0:
            #print('segment added')
            NRanges=np.concatenate((NRanges,Indexes[i:i+2]) ,axis=1)
    NRanges=np.delete(NRanges,0,1)
    return NRanges

def InitialComponents(a1,e1,a2,e2,Mstar,m1,m2,N):
    L = np.linspace(0, 2 * pi, N)   #testing
    R= np.zeros((2,N))  #radial collision points
    C = np.zeros((2, N))   #theta collision points

    for i in range(1,N):
        [[C[0,i], C[1,i]], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])


    #a data stuff
    ThetaKepler=np.sqrt((G*Mstar)/(np.power((R[:,:]),3)))
    Rdot1=rdot(a1,e1,C[:,:],m1+Mstar)
    ThetadotK1=thetadot(a1,e1,C[:,:],m1+Mstar)-ThetaKepler
    Rdot2 = rdot(a2, e2, C[:, :]-L, m2 + Mstar)
    ThetadotK2 = thetadot(a2, e2, C[:, :]-L, m2 + Mstar)-ThetaKepler
    abtitle=['a','b']
    RDotMultiple=np.zeros((2,N-1)) #Setting Up RDotMultiple
    ThetaDotMultiple = np.zeros((2, N - 1))  # Setting Up ThetaDotMultiple
    for i in range(0,2): #loop plotting for both collision points
    #rdot
        plt.figure(1+i)
        plt.subplot(221)
        plt.plot(L[1:N]/pi, Rdot1[i,1:N], label="Rdot1") #rdot in first orbit
        plt.plot(L[1:N]/pi, Rdot2[i,1:N], label="Rdot2") #rdot in the second orbit
        plt.legend() #attaching a legend to the graph
        plt.title('Rdot against Lambda/pi')
        plt.xlabel('Lambda/pi, Sepetation of the Orbits')
        plt.ylabel('Rdot')
        plt.xlim([0,2])
        plt.subplot(222) #plotting the rdots multiplied to tell the sign
        RDotMultiple[i,:]=np.multiply(Rdot1[i,1:N],Rdot2[i,1:N])
        plt.plot(L[1:N]/pi, RDotMultiple[i,:], label="Rdot1*Rdot2")
        plt.xlim([0,2])
        LRangeRDot = np.argwhere(np.diff(np.sign(RDotMultiple[i,:])) != 0) #Note these are still indices
        plotline = np.zeros((np.size(LRangeRDot), 1))
        plt.plot(((2) / N) * LRangeRDot, plotline, 'ro')

        NRangesRDotB = Segments(LRangeRDot, RDotMultiple[i, :])
        if i==0:
            NRangesRDotA=NRangesRDotB

        #athetadot
        plt.subplot(223)
        plt.plot(L[1:N]/pi, ThetadotK1[i,1:N], label="Thetadotk1")
        plt.plot(L[1:N]/pi, ThetadotK2[i,1:N], label="Thetadotl2")
        plt.legend()  # attaching a legend to the graph
        plt.title('Thetadotk against Lambda/pi')
        plt.xlabel('Thetadot')
        plt.xlim([0, 2 ])
        plt.subplot(224)
        ThetaDotMultiple[i,:]=np.multiply(ThetadotK1[i,1:N], ThetadotK2[i,1:N])
        plt.plot(L[1:N]/pi, ThetaDotMultiple[i,:], label="Rdot1*Rdot2")
        plt.xlim([0,2])
        LRangeThetaDot =  np.argwhere(np.diff(np.sign(ThetaDotMultiple[i,:])) != 0)
        plotline=np.zeros((np.size(LRangeThetaDot),1))
        plt.plot(((2) / N) *LRangeThetaDot, plotline, 'ro')

        NRangesThetaDotB = Segments(LRangeThetaDot, ThetaDotMultiple[i, :])
        if i == 0:
            NRangesThetaDotA = NRangesThetaDotB

        plt.suptitle('collision point %s'%(abtitle[i]))





    #graph of the collision points

    plt.figure(3)
    plt.plot(L[1:N] / pi, R[0, 1:N], label="Collision RAdius a")
    plt.plot(L[1:N] / pi, R[1, 1:N], label="Collision Radius b")
    plt.plot([0,2],[1,1], label="1au")
    plt.legend()
    plt.title('Collision Points')

   #finding min values from less than 1au

    LRangeDistanceA= np.argwhere(np.diff(np.sign(R[0,1:N] - 1)) != 0)  #the points at which the graph changes sign
    LRangeDistanceB  = np.argwhere(np.diff(np.sign(R[1, 1:N] - 1)) != 0)
    plt.plot(((2)/N)*LRangeDistanceA,[1,1],'ro')
    plt.plot(((2) / N) * LRangeDistanceB, [1, 1], 'bo')
    plt.show(3)

     #creating segments

    #Orbit A
    NRangesDistanceA=Segments(LRangeDistanceA, (R[0,:]-1))
    print('NRangesDistanceA'), print(((2) / N) *NRangesDistanceA)
    print('NRagnesRDotA'), print(((2) / N) * NRangesRDotA)
    print('NRangesThetaDotA'), print(((2) / N) *NRangesThetaDotA)

    #Orbit B
    NRangesDistanceB = Segments(LRangeDistanceB, (R[1, :] - 1))
    print('NRangesDistanceB'), print(((2) / N) * NRangesDistanceB)
    print('NRagnesRDotB'), print(((2) / N) * NRangesRDotB)
    print('NRangesThetaDotB'), print(((2) / N) * NRangesThetaDotB)
    #WANT TO RETURN AN OUTPUT OF a range of L for each [,] , [,]
    return

InitialComponents(2*au,0.99,2.1*au,0.993,1.2e30,2e10,2e10,1000)


