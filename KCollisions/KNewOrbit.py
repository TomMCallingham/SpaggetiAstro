from Funcs import *
from math import *
import matplotlib.pyplot as plt
G =6.674e-11# 2.982e-27
au= 1.496e8

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
'''
def NewOrbit(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,c1,r1): #This function gives the new orbit created, given two orbits and there crossing points
    m3=(m1+m2)   #for now assuming its perfectly merging
    #Calculating Constants
    u1=G*(m1+Mstar)
    u2=G*(m2+Mstar)
    u3=G*(m3+Mstar)
    A1=A(a1,e1,m1+Mstar)
    A2=A(a2,e2,m2+Mstar)
    A3=(u3*m3)/(((u1*m1)/A1)+((u2*m2)/A2))
    alpha=(u3/((A3**2.)*r1)) -1
    beta=((m1*rdot(a1,e1,c1-s1,Mstar+m1))+(m2*rdot(a2,e2,c1-s2,Mstar+m2)))/m3  #RDOT IS WRONG!
    #Calculating the New Orbit
    e3=sqrt((alpha**2.)+((beta**2.)/(A3**2.)))
    #if e3>=1: #Checking the new orbit is bound, giving an error if not
        #print('ERROR: Unbound Orbit!')
    a3=u3/((1-(e3**2.))*(A3**2.))
    s3= -acos(alpha/e3)  +c1 #note changed
    if abs(sin(c1-s3)-beta/(A3*e3))>10e-10: #checking the degeneracy of arccos
        s3=2*pi-(s3-2*c1)
    if abs(sin(c1-s3)-beta/(A3*e3))>10e-10: #checking the degeneracy of arccos
        print('Error in new orbit')
    #Outputting the NewOrbit Data
    NewOrbitData=np.array([a3,e3,s3,m3])
    return NewOrbitData
'''
'''
def KNewOrbit(a1,e1,s1,a2,e2,s2,Mstar,m1,m2,K,c1,r1): #This function gives the new orbit created, given two orbits and there crossing points
    m3=K*(m1+m2)   #for now assuming its perfectly merging
    #Calculating Constants
    u1=G*(m1+Mstar)
    u2=G*(m2+Mstar)
    u3=G*(m3+Mstar)
    A1=A(a1,e1,m1+Mstar)
    A2=A(a2,e2,m2+Mstar)
    A3=(u3*m3)/(((u1*m1)/A1)+((u2*m2)/A2))
    alpha=(u3/((A3**2.)*r1)) -1
    beta=((m1*rdot(a1,e1,c1-s1,Mstar+m1))+(m2*rdot(a2,e2,c1-s2,Mstar+m2)))/m3  #RDOT IS WRONG!
    #Calculating the New Orbit
    e3=sqrt((alpha**2.)+((beta**2.)/(A3**2.)))
    #if e3>=1: #Checking the new orbit is bound, giving an error if not
        #print('ERROR: Unbound Orbit!')
    a3=u3/((1-(e3**2.))*(A3**2.))
    s3= -acos(alpha/e3)  +c1 #note changed
    if abs(sin(c1-s3)-beta/(A3*e3))>10e-10: #checking the degeneracy of arccos
        s3=2*pi-(s3-2*c1)
    if abs(sin(c1-s3)-beta/(A3*e3))>10e-10: #checking the degeneracy of arccos
        print('Error in new orbit')
    #Outputting the NewOrbit Data
    NewOrbitData=np.array([a3,e3,s3,m3])
    return NewOrbitData
'''

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
    plt.plot(L[1:1000],R[0,1:1000],label="Collision Point 1")
    plt.plot(L[1:1000], R[1, 1:1000], label="Collision Point 2")
    plt.plot([0,2*pi],[1,1], label="1 au")
    plt.legend() #attaching a legend to the graph
    plt.title('Radius of the Collision Points against Lambda')
    plt.xlabel('Lambda, Sepetation of the Orbits')
    plt.ylabel('Radius of the the Collision Points (au)')
    plt.xlim([0,2*pi])
    plt.show()
    return





