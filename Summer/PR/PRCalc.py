import matplotlib.pyplot as plt
import numpy as np
from Functions import *
Lsun=3.828e26
#Values
a1 = 2.5 * au
e1 = 0.998
a2 = 2.4 * au
e2 = 0.997
I=1* ((2 * pi) / 360)

def Luminosity(t): #note time in years
    L=3.26*Lsun*((0.1+t*(10**-6.))**(-1.18))
    return L
def GraphLuminosity():
    Time=np.linspace(0,9*(10**9.),1e3)
    L=Luminosity(Time)
    plt.figure()
    plt.loglog(Time,L/Lsun,label='Luminosity in Lsun')
    plt.ylabel('Luminosity (Lsun)')
    plt.xlabel('Time')
    plt.title('WD Luminosity against time')
    plt.legend()
    plt.show()
    return

