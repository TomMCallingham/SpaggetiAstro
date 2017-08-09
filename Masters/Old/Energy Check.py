from Funcs import *
import numpy as np
import matplotlib.pyplot as plt

def EnergyCheck(a,e,Mstar):
    F = np.linspace(0, 2 * pi, 1000)   #NOTE CHECK L def!! + or -
    R=npr(a,e,F)
    Rdot=rdot(a,e,F,Mstar)
    Thetadot=thetadot(a,e,F,Mstar)
    E=SEnergy(Rdot, Thetadot, R, Mstar)
    EAlt = -(G * Mstar) / (2 * a)
    plt.figure()
    plt.plot(F/pi,E, label="EnergyFunc")
    plt.plot([0,2],[EAlt, EAlt], label="Alternate Energy")
    plt.legend()
    plt.show()
    return

EnergyCheck(2 , 0.99 , 1.2e30)