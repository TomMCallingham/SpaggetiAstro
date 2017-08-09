from Funcs import *
from math import *
import matplotlib.pyplot as plt
G =6.674e-11# 2.982e-27
au= 1.496e8
Mstar=1.2e30

def CollisionPoints(a1,e1,s1,a2,e2,s2): #This function finds the crossing points of two given ellipse,
    #function will give an error if L=0,e2=0
    #Calculating Constants used in the method
    L=s1-s2
    k=(a2*(1-(e2**2)))/(a1*(1-(e1**2)))
    alpha=(cos(L)-k*(e1/e2))/sin(L) #note e2=/=0
    beta=(1-k)/(e2*sin(L))
    #Calculating the two soluions:
    #the first solution
    c1= acos((-(alpha * beta) + sqrt(((alpha * beta) ** 2) - (((alpha ** 2) + 1) * ((beta ** 2) - 1)))) / ((alpha ** 2) + 1)) + s1
    if abs(alpha*cos(c1-s1)+beta-sin(c1-s1))>10e-6 : #checking the degeneracy of arccos
        c1 = 2*pi -(c1-2*s1)
    r1=npr(a1,e1,c1-s1)
    #the second solution
    c2= acos((-(alpha * beta) - sqrt(((alpha * beta) ** 2) - (((alpha ** 2) + 1) * ((beta ** 2) - 1)))) / ((alpha ** 2) + 1)) + s1
    if abs(alpha*cos(c2-s1)+beta-sin(c2-s1))>10e-6:
        c2 = 2*pi - (c2-2*s1)
    r2=npr(a1,e1,c2-s1)
    #Output is the two collision points in angle and radius in an array
    CollisionData=np.array([[c1,c2],[r1,r2]])
    return CollisionData

def NewOrbit(a1, e1, s1, a2, e2, s2, Mstar, m1, m2, C, R): #This function gives the new orbit created, given two orbits and there crossing points
    m3=(m1+m2)   #for now assuming its perfectly merging
    #Calculating Constants
    rd1=rdot(a1, e1, C - s1, Mstar)
    td1=thetadot(a1, e1, C - s1, Mstar)
    rd2=rdot(a2, e2, C - s2, Mstar)
    td2=thetadot(a2, e2, C - s2, Mstar)
    rd3=((m1*rd1)+(m2*rd2))/m3
    td3 = ((m1 * td1) + (m2 * td2)) / m3
    e3=newe(rd3,td3,R,Mstar)
    a3=newa(rd3,td3,R,Mstar)

    s3=C- np.arccos((((a3/R)*(1-(e3**2.)))-1)/e3)
    '''
    print('e3',e3)
    print('a3', a3)
    print('s3', s3)
    '''
    if abs(rd3-((e3*((G*Mstar))*sin(C-s3))/((R**2.)*td3)))>10e-6: #checking the degeneracy of arccos
        s3=2*C-2*pi-s3
    if abs(rd3-((e3*((G*Mstar))*sin(C-s3))/((R**2.)*td3)))>10e-6: #checking the degeneracy of arccos
        print('Error in new orbit')
    '''print('s3fix', s3)'''
    #Outputting the NewOrbit Data
    NewOrbitData=np.array([a3,e3,s3,m3])
    return NewOrbitData

def CollisionGraph(a1,e1,s1,a2,e2,s2,Mstar,m1,m2):
    #Setting up the Plotter
    F = np.linspace(0, 2 * pi, 100)
    plt.figure()
    #Plotting 1
    R = npr(a1, e1, F - s1)
    plt.polar(F, R,label="Parent Orbit 1")
    #Plotting 2
    R = npr(a2, e2, F - s2)
    plt.polar(F, R, label="Parent Orbit 2")
    #Plotting the first New Orbit
    CollisionData=CollisionPoints(a1,e1,s1,a2,e2,s2) #Loading the ellipse intersection points
    plt.scatter(CollisionData[0,:], CollisionData[1,:], c='r') #Highlighting the collision points
    [a3, e3, s3,m3]=NewOrbit(a1, e1, s1, a2, e2, s2,Mstar,m1, m2,CollisionData[0,0],CollisionData[1,0]) #loading the data on the new orbit
    if e3<1 : #Checking its a valid bound orbit
        R = npr(a3, e3, F - s3)
        plt.polar(F, R, label="New Orbit 1")
    else: #else give an error
        print('Error: Orbit 3 Unbound!')
    #Plotting the Second New Orbit
    [a4, e4, s4,m4] = NewOrbit(a1, e1, s1, a2, e2, s2, Mstar, m1, m2, CollisionData[0,1], CollisionData[1, 1])
    if e4 < 1:
        R = npr(a4, e4, F - s4)
        plt.polar(F, R, label="New Orbit 2")
    else:
        print('Error: Orbit 4 Unbound!')
    plt.legend() #attaching a legend to the graph
    #plt.title('Graphing the Bound Orbits with K=%s'%K)
    plt.show()
    return

def RadialGraph(a1,e1,a2,e2):
    s1=0
    L = np.linspace(0, 2 * pi, 1000)
    R= np.zeros((2,1000))
    for i in range(1,1000):
        [[c1, c2], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, -L[i])
     #now in au Units, no scaling needed (1/au)*R
    '''
    plt.plot(L[1:1000],R[0,1:1000]/au,label="Collision Point a")
    plt.plot(L[1:1000], R[1, 1:1000]/au, label="Collision Point b")
    '''



    plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    plt.rcParams.update({'font.size': 12})
    plt.semilogy(L[1:1000]/pi, R[0, 1:1000] / au, label="Collision Point A")
    plt.semilogy(L[1:1000]/pi, R[1, 1:1000] / au, label="Collision Point B")
    plt.semilogy([0,2],[0.013,0.013], label='Maximum Roche Radius')
    plt.semilogy([0, 2], [0.005, 0.005], label='Solar Radius')
    #plt.plot([0,2*pi],[1,1], label="1 au")
    plt.legend() #attaching a legend to the graph
    #plt.title('Radius of the Collision Points against Lambda')
    plt.title(r"$r_ {*}$" + ' against ' + r"$\lambda$")

    plt.xlabel(r"$\lambda$" +'/pi' )
    plt.ylabel(r"$r_ {*}$"+' /au')
    plt.xlim([0,2])
    plt.show()
    return

def SequenceCollisionGraph(ap1,ep1,sp1,ap2,ep2,sp2):
    x = 'a'
    CollisionData = CollisionPoints(ap1, ep1, sp1, ap2, ep2, sp2)
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    shift = (-pi / 2) - C
    sp1 = sp1 + shift
    sp2 = sp2 + shift
    plt.figure(1)




    #Setting up the Plotter
    m1=1
    m2=1
    F = np.linspace(0, 2 * pi, 100)

    #Plotting 1
    R1= npr(ap1, ep1, F - sp1)
    plt.polar(F, R1/au,label="Parent Ring 1", c='darkorange')
    #Plotting 2
    R2 = npr(ap2, ep2, F - sp2)
    plt.polar(F, R2/au, label="Parent Ring 2",c='g')
    #Plotting the first New Orbit
    CollisionData=CollisionPoints(ap1,ep1,sp1,ap2,ep2,sp2) #Loading the ellipse intersection points

    plt.scatter(0,0,c='brown', marker='*', zorder=9)
    plt.scatter(CollisionData[0, 0], CollisionData[1, 0] / au, c='r', zorder=10)  # Highlighting the collision points
    #plt.scatter(CollisionData[0, 1], CollisionData[1, 1] / au, c='r', zorder=10)  # Highlighting the collision points


    plt.polar([pi+sp2,sp2],[np.max(R2)/au,np.max(R2)/(2*au)],c='k')
    plt.polar([pi + sp1, sp1], [np.max(R1) / au, np.max(R1) / (2 * au)], c='k')
    plt.polar(np.linspace(sp1,sp2,100), 1*np.ones(100), c='m',label='Lambda')
    plt.polar([pi/2,-pi/2],[3,3],'b', label='Intersection Line')

    
    #plt.scatter(CollisionData[0, 1], CollisionData[1, 1], c='r')  # Highlighting the collision points

    #First Parent Graph
    '''
    #Gen1
    [ac1, ec1, sc1,m3]=NewOrbit(ap1, ep1, sp1, ap2, ep2, sp2,Mstar,m1, m2,CollisionData[0,0],CollisionData[1,0]) #loading the data on the new orbit
    R = npr(ac1, ec1, F - sc1)
    plt.polar(F, R/au, label="Gen1:C1=P1+P2")


    #Gen 2
    [ac2, ec2, sc2, m4] = NewOrbit(ap1, ep1, sp1, ac1, ec1, sc1, Mstar, m1, m2, CollisionData[0, 0],
                                CollisionData[1, 0])  # loading the data on the new orbit
    R = npr(ac2, ec2, F - sc2)
    plt.polar(F, R/au, label="Gen2:C2=P1+C1")
    [ac3, ec3, sc3, m4] = NewOrbit(ap2, ep2, sp2, ac1, ec1, sc1, Mstar, m1, m2, CollisionData[0, 0],
                                CollisionData[1, 0])  # loading the data on the new orbit
    R = npr(ac3, ec3, F - sc3)
    plt.polar(F, R/au, label="Gen2:C3=P2+C1")

    # Gen 3
    [ac4, ec4, sc4, m4] = NewOrbit(ap1, ep1, sp1, ac2, ec2, sc2, Mstar, m1, m2, CollisionData[0, 0],
                                CollisionData[1, 0])  # loading the data on the new orbit
    R = npr(ac4, ec4, F - sc4)
    plt.polar(F, R/au, label="Gen3:C4=P1+C2")

    [ac5, ec5, sc5, m4] = NewOrbit(ap1, ep1, sp1, ac2, ec2, sc2, Mstar, m1, m2, CollisionData[0, 0],
                                CollisionData[1, 0])  # loading the data on the new orbit
    R = npr(ac5, ec5, F - sc5)
    plt.polar(F, R/au, label="Gen3:C5=P1+C3")

    [ac6, ec6, sc6, m4] = NewOrbit(ac1, ec1, sc1, ac2, ec2, sc2, Mstar, m1, m2, CollisionData[0, 0],
                                CollisionData[1, 0])  # loading the data on the new orbit
    R = npr(ac6, ec6, F - sc6)
    plt.polar(F, R/au, label="Gen3:C6=C1+C2")

    [ac7, ec7, sc7, m4] = NewOrbit(ac1, ec1, sc1, ac3, ec3, sc3, Mstar, m1, m2, CollisionData[0, 0],
                                CollisionData[1, 0])  # loading the data on the new orbit
    R = npr(ac7, ec7, F - sc7)
    plt.polar(F, R/au, label="Gen3:C7=C1+C3")

    [ac8, ec8, sc8, m4] = NewOrbit(ac2, ec2, sc2, ac3, ec3, sc3, Mstar, m1, m2, CollisionData[0, 0],
                                CollisionData[1, 0])  # loading the data on the new orbit
    R = npr(ac8, ec8, F - sc8)
    plt.polar(F, R/au, label="Gen3:C8=C2+C3")



    '''
    plt.title('Orbits in au')
    plt.legend(loc=2)
    plt.show()
    #print([ac1/au, ec1, sc1])
    return
SequenceCollisionGraph(2*au,0.9,0,2.1*au,0.8,1)
'''
RadialGraph(2.5*au,0.998,2.4*au,0.997)

Rwd=3e-6
q=200*Rwd
emin=1-(q/2.5)
print('emin',emin)
'''
#RadialGraph(2.5*au,0.999,2.5*au,0.9998)
