from math import *
import numpy as np
G=2.982e-27  #now in au, before 6.67408e-11 #Graviational Constant


def rdot(a,e,f,m):
    return np.sqrt((G*m)/(a*(1-(e**2))))*e*np.sin(f)

def thetadot(a,e,f,m):
    return np.sqrt((G*m*a*(1-(e**2.))))/(npr(a,e,f)**2)


def npr(a,e,f): #Numpyr, used for plotting
    return (a*(1-(e**2)))/(1+(e*np.cos(f)))

def newe(Rdot,Thetadot,R,Mstar):
    u=G*(Mstar)  #could have Mstar+m3, but not importatnt
    e=np.sqrt((((((R**3.)*(Thetadot**2.))/u)-1)**2.)+(((Rdot**2.)*(R**4)*(Thetadot**2))/(u**2.)))
    return e

def SEnergy(rdot, thetadot, r, Mstar):
    E=((((r*thetadot)**2.)+(rdot**2.))/2)-((G*Mstar)/r)
    return E

def thetakepler(R,Mstar):
    ThetaK = np.sqrt((G * Mstar) / (R ** 3.))
    return ThetaK


def hangle(a,e,m):
    return np.sqrt(G*m*a*(1-(e**2.)))
#def thetadot(a,e,f,m):
  #  return hangle(a,e,m)/(npr(a,e,f)**2.)
#def rdot(a, e, f, m):
   # return A(a, e, m) * e * np.sin(f)

def A(a, e, m):
    Acalc = np.sqrt((G * m) / (a * (1 - (e ** 2.))))
    return Acalc


def r(a, e, f):  # note f is from the pericentere, not reference
    return (a * (1 - (e ** 2))) / (1 + (e * cos(f)))

