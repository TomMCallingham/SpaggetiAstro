from PRDrag.EvolvePRe import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def ringPR(a,e,Tini):
    Tyears=10**8.
    TimeStep=100
    plt.figure()

    F = np.linspace(0, 2 * pi, 1000)
    # Plotting 1
    colours = cm.rainbow(np.linspace(0, 1, 5))
    R1 = npr(a, e, F)
    plt.polar(F, R1 / au, label='initial')
    for i in range(-2,2):
        Rfrag=10**i
        (Data, Tmini)=PREvolve(a,e,Tyears,Tini,Rfrag)
        T=np.size(Data[0,:])-1
        for t in range(0,T+1):
            RelativeAngle='isthisevenneeded'

        for t in (0,4): #range(0,5)
            Timeshow=int(T*(10**-t))
            R1 = npr(Data[1,Timeshow], Data[2,Timeshow], F)
            plt.polar(F,R1/au,label='Rfrag=10^%s, T=10^%s'%(i,6-t), color=colours[t])



    plt.legend()

   # plt.loglog([1,Tyears],[10**3.,10**3.], label='Km')


    plt.show()
    print('finished graphing')
    return

#changeperorbit(5e9,2.5*au,0.995)
#ringdispresion(2.1*au,0.995,10**7.,0.5e9)
ringPR(2.5 * au, 0.995, 0.5e9)
#GraphPrEvolve(2.5*au,0.995,10**8,0.5e9)

'''
        Angle=np.zeros((1,T+1))
        for t in range(0,T+1):
            Angle
        '''