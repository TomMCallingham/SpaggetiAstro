from math import *
import numpy as np
G=6.67408e-11 #Graviational Constant
def A(a, e, m):
    Acalc=sqrt((G*m)/(a*(1-(e**2.))))
    return Acalc
def r(a,e,f): #note f is from the pericentere, not reference
    return (a*(1-(e**2)))/(1+(e*cos(f)))
def rdot(a,e,f,m):
    return A(a,e,m)*e*sin(f)
def hangle(a,e,m):
    return sqrt(G*m*a*(1-(e**2.)))
def thetadot(a,e,f,m):
    return hangle(a,e,m)/(r(a,e,f)**2.)
#Numpyr, used for plotting
def npr(a,e,f): #note f is from the pericentere, not reference
    return (a*(1-(e**2)))/(1+(e*np.cos(f)))





