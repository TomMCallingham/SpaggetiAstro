import numpy as np
import matplotlib.pyplot as plt
def vEnBoundary(Ehat,):
    return


def dKFunc(v1,m1,v2,m2):
    zmf=((m1*v1)+(m2*v2))/(m1+m2)
    Ebefore=((m1*(v1**2.))+(m2*(v2**2.)))/2
    VColl=abs(v1-v2)
    SpecE=((VColl**2.)*(m2/m1))/2   #CHECK

    Kn=10
    K=np.linspace(0,1,Kn+1)
    K=K[1:Kn+1]
    Vmax=np.sqrt(((2*Ebefore)/(m1+m2))*((1-K)/K))

    plt.plot(Vmax,K)
    plt.xlabel('Vmax')
    plt.ylabel('K')
    plt.show()
    return


dKFunc(10,1,20,1)