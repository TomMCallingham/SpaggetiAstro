from math import *
import numpy as np
G =6.674e-11# 2.982e-27
au= 1.496e11
year = 3.154e+7
Mstar = 1.2e30

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

def newa(rd,td,R,Mstar):
    e=newe(rd,td,R,Mstar)
    u = G * (Mstar)
    a=((td**2.)*(R**4.))/(u*(1-(e**2.)))
    return a


def SEnergy(rdot, thetadot, r, Mstar):
    E=((((r*thetadot)**2.)+(rdot**2.))/2)-((G*Mstar)/r)
    return E

def thetakepler(R,Mstar):
    ThetaK = np.sqrt((G * Mstar) / (R ** 3.))
    return ThetaK


def hangle(a,e,m):
    return np.sqrt(G*m*a*(1-(e**2.)))


def Timeperiod(a): #a in au
    a=a*au
    T=(2*pi)*np.sqrt((a**3.)/(G*Mstar))
    T=T/year
    #print(T)
    return
#def thetadot(a,e,f,m):
  #  return hangle(a,e,m)/(npr(a,e,f)**2.)
#def rdot(a, e, f, m):
   # return A(a, e, m) * e * np.sin(f)
'''
def A(a, e, m):
    Acalc = np.sqrt((G * m) / (a * (1 - (e ** 2.))))
    return Acalc


def r(a, e, f):  # note f is from the pericentere, not reference
    return (a * (1 - (e ** 2))) / (1 + (e * cos(f)))

'''
#Timeperiod(2.5)
#Timeperiod(2.4)