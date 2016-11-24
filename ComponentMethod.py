import numpy as np
from math import *
from Funcs import *
from NewOrbit import *

au= 149597871e3
G = 6.67408e-11
#big data table of N1 and N2 values, storing eccentriciy
def newmom(dot1,dot2,m1,m2,N1,N2):
    T=((N1*m1*dot1)+(N2*m2*dot2))/((N1*m1)+(N2*m2))
    return T
def newe(Rdot,Thetadot,R,m3,Mstar):
    u=G*(Mstar+m3)
    e=sqrt((((((R**3.)*(Thetadot**2.))/u)-1)**2.)+(((Rdot**2.)*(R**4)*(Thetadot**2))/(u**2.)))
    return e
def CascadeComponents(a1,e1,s1,a2,e2,s2,Mstar,m1,m2):
    N=3
    [[ca, cb], [ra, rb]] = CollisionPoints(a1, e1, s1, a2, e2, s2)
    ra1dot = rdot(a1, e1, ca - s1, m1)  # CHECK SIGN OF ANGLES
    thetaa1dot=thetadot(a1, e1, ca - s1, m1)
    ra2dot = rdot(a2, e2, ca - s2, m2)
    thetaa2dot = thetadot(a2, e2, ca - s2, m2)
    rb1dot = rdot(a1, e1, cb - s1, m1)
    thetab1dot = thetadot(a1, e1, cb - s1, m1)
    rb2dot = rdot(a2, e2, cb - s2, m2)
    thetaa1dot = thetadot(a2, e2, cb - s2, m2)
    #print('Orbit_a'),print('1'), print(ra1dot), print(thetaa1dot), print('2'), print(ra2dot), print(thetaa2dot)
    #print('Orbit_b'), print(rb1dot), print(rb2dot)

    #a for now
    aData=np.zeros((N+1,N+1,4))
    #ParentOrbit 1
    aData[1, 0, 1] =  m1  # new m
    aData[1, 0, 2] = ra1dot
    aData[1, 0, 3] = thetaa1dot
    aData[1, 0, 0] = e1
    #Parent Orbit 2
    aData[0, 1, 1] =  m2  # new m
    aData[0, 1, 2] = ra2dot
    aData[0, 1, 3] = thetaa2dot
    aData[0, 1, 0] = e2
    for N1 in range(1,N+1):
        for N2 in range(1,N+1):
            aData[N1,N2,1]=N1*m1+N2*m2 #new m
            aData[N1,N2,2]=newmom(ra1dot,ra2dot,m1,m2,N1,N2) #new rdot
            aData[N1,N2,3]=newmom(thetaa1dot,thetaa2dot,m1,m2,N1,N2) #new thetadot
            aData[N1,N2,0]=newe(aData[N1,N2,2],aData[N1,N2,3],ra,aData[N1,N2,1],Mstar) #new e
    return aData

Data=CascadeComponents(2*au,0.99,0,2.1*au,0.993,2,1.2e30,2e10,2e10)
print(Data[:,:,0])