from math import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as npal


G =6.674e-11# 2.982e-27
au= 1.496e11
year = 3.154e+7
Mstar = 1.2e30
density=2e3 #Kg


def npr(a,e,f): #Numpyr, used for plotting
    return (a*(1-(e**2)))/(1+(e*np.cos(f)))

def rdot(a,e,f): #radial velocity
    return np.sqrt((G*Mstar)/(a*(1-(e**2))))*e*np.sin(f)

def thetadot(a,e,f): #angular velocity
    return np.sqrt((G*Mstar*a*(1-(e**2.))))/(npr(a,e,f)**2)




def CollisionPoints(a1,e1,s1,a2,e2,s2): #This function finds the crossing points of two given ellipse,
    #function will give an error if L=0,e2=0
    #Calculating Constants used in the method
    L=s1-s2 #IMPORTANT DEF
    k=(a2*(1-(e2**2)))/(a1*(1-(e1**2)))
    alpha=(cos(L)-k*(e1/e2))/sin(L) #note e2=/=0
    beta=(1-k)/(e2*sin(L))
    #Calculating the two soluions:
    #the first solution
    delta=((alpha * beta) ** 2) - (((alpha ** 2) + 1) * ((beta ** 2) - 1))
    if delta>=0:
        cosc1 = (-(alpha * beta) + sqrt(delta)) / (
        (alpha ** 2) + 1)
        if abs(cosc1)<=1:
            c1= acos(cosc1) + s1
            if abs(alpha*cos(c1-s1)+beta-sin(c1-s1))>10e-6 : #checking the degeneracy of arccos
                c1 = 2*pi -(c1-2*s1)
            r1=npr(a1,e1,c1-s1)
        else:
            #print('No A Point')
            c1 = 0
            r1 = 0

    #the second solution
        cosc2=(-(alpha * beta) - sqrt(delta)) / ((alpha ** 2) + 1)
        if abs(cosc2)<=1:

            c2= acos(cosc2) + s1
            if abs(alpha*cos(c2-s1)+beta-sin(c2-s1))>10e-6:
                c2 = 2*pi - (c2-2*s1)
            r2=npr(a1,e1,c2-s1)
        else:
            #print('No B point')
            c2=0
            r2=0
    else:
        #print('No points')
        #print('L=',L)
        c1 = 0
        r1 = 0
        c2=0
        r2=0
    #Output is the two collision points in angle and radius in an array
    CollisionData=np.array([[(c1)% (2 * pi),(c2)% (2 * pi)],[r1,r2]])  #r1 is a or 0, r1 is b or 1. r2>r1
    return CollisionData



def newe(Rdot,Thetadot,R):
    u=G*(Mstar)
    e=np.sqrt((((((R**3.)*(Thetadot**2.))/u)-1)**2.)+(((Rdot**2.)*(R**4)*(Thetadot**2))/(u**2.)))
    return e

def newa(rd,td,R):
    e=newe(rd,td,R)
    u = G * (Mstar)
    a=((td**2.)*(R**4.))/(u*(1-(e**2.)))
    return a

def orbitalvalues(rsource):
    p=0.01
    a=rsource/2
    e=1-p/a
    #print('rsource=',rsource,a,e )
    return(a*au,e)

def orbitalvaluesmod(rsource,pmod):
    p=0.01*pmod
    a=rsource/2
    e=1-p/a
    #print('rsource=',rsource,a,e )
    return(a*au,e)