from Funcs import *
from math import *
import matplotlib.pyplot as plt
G = 2.982e-27#6.67408e-11
au= 1#149597871e3

def CollisionPoints(a1,e1,s1,a2,e2,s2): #This function finds the crossing points of two given ellipse
    #function will give an error if L=0,e2=0
    #Calculating Constants used in the method
    L=s1-s2
    k=(a2*(1-(e2**2)))/(a1*(1-(e1**2)))
    alpha=(cos(L)-k*(e1/e2))/sin(L) #note e2=/=0
    beta=(1-k)/(e2*sin(L))
    #Calculating the two soluions:
    #the first solution
    c1= acos((-(alpha * beta) + sqrt(((alpha * beta) ** 2) - (((alpha ** 2) + 1) * ((beta ** 2) - 1)))) / ((alpha ** 2) + 1)) + s1
    if abs(alpha*cos(c1-s1)+beta-sin(c1-s1))>10e-10 : #checking the degeneracy of arccos
        c1 = 2*pi -(c1-2*s1)
    r1=r(a1,e1,c1-s1)
    #the second solution
    c2= acos((-(alpha * beta) - sqrt(((alpha * beta) ** 2) - (((alpha ** 2) + 1) * ((beta ** 2) - 1)))) / ((alpha ** 2) + 1)) + s1
    if abs(alpha*cos(c2-s1)+beta-sin(c2-s1))>10e-10:
        c2 = 2*pi - (c2-2*s1)
    r2=r(a1,e1,c2-s1)
    #Output is the two collision points in angle and radius in an array
    CollisionData=np.array([[c1,c2],[r1,r2]])
    return CollisionData

def EqualNewOrbit(a1,e1,s1,a2,e2,s2,Mstar,c1,r1): #This function gives the new orbit created, given two orbits and there crossing points
    #Calculating Constants
    u1=G*(Mstar)
    u2=G*(Mstar)
    u3=G*(Mstar)
    A1=A(a1,e1,Mstar)
    A2=A(a2,e2,Mstar)
    A3=(u3*2)/((u1/A1)+(u2/A2)) #NOO
    alpha=(u3/((A3**2.)*r1)) -1
    beta=(rdot(a1,e1,c1-s1,Mstar)+rdot(a2,e2,c1-s2,Mstar))/2  #RDOT
    #Calculating the New Orbit
    e3=sqrt((alpha**2.)+((beta**2.)/(A3**2.)))
    if e3>=1: #Checking the new orbit is bound, giving an error if not
        print('ERROR: Unbound Orbit!')
    a3=u3/((1-(e3**2.))*(A3**2.))
    s3= acos(alpha/e3)+c1
    if abs(sin(s3+c1)-beta/(A3*e3))>10e-10: #checking the degeneracy of arccos
        s3=2*pi-(s3-2*c1)
    #Outputting the NewOrbit Data
    NewOrbitData=np.array([a3,e3,s3])
    return NewOrbitData

def EqualCollisionGraph(a1,e1,s1,a2,e2,s2,Mstar):
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
    [a3, e3, s3]=EqualNewOrbit(a1, e1, s1, a2, e2, s2,Mstar,CollisionData[0,0],CollisionData[1,0]) #loading the data on the new orbit
    if e3<1 : #Checking its a valid bound orbit
        R = npr(a3, e3, F - s3)
        plt.polar(F, R, label="New Orbit 1")
    else: #else give an error
        print('Error: Orbit 3 Unbound!')
    #Plotting the Second New Orbit
    [a4, e4, s4] = EqualNewOrbit(a1, e1, s1, a2, e2, s2, Mstar, CollisionData[0,1], CollisionData[1, 1])
    if e4 < 1:
        R = npr(a4, e4, F - s4)
        plt.polar(F, R, label="New Orbit 2")
    else:
        print('Error: Orbit 4 Unbound!')
    plt.legend() #attaching a legend to the graph
    plt.title('Graphing the Bound Orbits')
    plt.show()
    return

def EqualRadialGraph(a1,e1,a2,e2):
    s1=0
    L = np.linspace(0, 2 * pi, 1000)
    R= np.zeros((2,1000))
    for i in range(1,1000):
        [[c1, c2], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, -L[i])
    R=(1/au)*R
    plt.plot(L[1:1000],R[0,1:1000],label="Collision Point 1")
    plt.plot(L[1:1000], R[1, 1:1000], label="Collision Point 2")
    plt.plot([0,2*pi],[1,1], label="1 au")
    plt.legend() #attaching a legend to the graph
    plt.title('Radius of the Collision Points against Lambda')
    plt.xlabel('Lambda, Sepetation of the Orbits')
    plt.ylabel('Radius of the the Collision Points (au)')
    plt.show()
    return




#EqualCollisionGraph(2*au,0.3,0,2.1*au,0.1,2,1.2e30)
#plt.show()
