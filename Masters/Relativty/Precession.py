import numpy as np
from Masters.Funcs import *
import matplotlib.pyplot as plt
Mstar = 1.2e30
au = 1.496e11
year = 3.154e+7
#Setup
x='a'
a1=2.5*au
e1=0.998
s1=0
a2=2.4*au
e2=0.997
s2=1
def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5) #in MYrs
    wp = (2 * np.pi) / Tp
    return (Tp,wp)

(Tp1,wp1)=TPrecess(a1,e1)
(Tp2,wp2)=TPrecess(a2,e2)


print('Tp1/MYrs',Tp1)
print('Tp2/Myrs',Tp2)
print('wp1*Myrs',wp1)
print('wp2*Myrs',wp2)

'''
Dwp=abs(wp1-wp2)
DTp=(2*np.pi)/Dwp
print('DTp',DTp)
print('Dwp',Dwp)
print('Tp2/Tp1',Tp2/Tp1)
'''

def intersectgraphs():
    s1i=0
    s2i=0
    T=np.arange(0,20,10**-4.) #megayears
    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)

    R1plus=npr(a1,e1,(pi/2)+s1i+wp1*T)
    R1minus = npr(a1, e1, -(pi / 2) + s1i + wp1 * T)
    R2plus = npr(a2, e2, (pi / 2) + s2i + wp2 * T)
    R2minus = npr(a2, e2, -(pi / 2) + s2i + wp2 * T)

    '''
    #r graph
    plt.figure()
    plt.subplot(121)
    plt.plot(T,R1plus/au,label='R1plus')
    plt.plot(T,R2plus/au,  label='R2plus')
    plt.ylabel('radius at plane intersect in au')
    plt.xlabel('Time in Mega Years')
    plt.legend()
    plt.title('plus')
    plt.subplot(122)
    plt.plot(T,R1minus/au,  label='R1minus')
    plt.plot(T,R2minus/au,  label='R2minus')
    #plt.xlabel('Time in Mega Years')
    #plt.ylabel('radius at plane intersect in au')
    plt.legend()
    plt.title('minus')
    plt.suptitle('(a1=%s,e1=%s), (a2=%s,e2=%s)'%(a1/au, e1, a2/au, e2))
    '''

    #Log
    plt.figure()
    #plt.subplot(211)

    plt.semilogy(T, R1plus / au, label='R1plus')
    plt.semilogy(T, R2plus / au, label='R2plus')
    plt.ylabel('Radius at plane intersect in au')
    plt.title('Radius of Orbits on the Line of Intersection')
    plt.xlabel('Time in Mega Years')
    plt.legend()
    #plt.title('plus')

   # plt.subplot(212)
    plt.figure()
    plt.semilogy(T, R1minus / au, label='R1minus')
    plt.semilogy(T, R2minus / au, label='R2minus')
    plt.xlabel('Time in Mega Years')
    #plt.xlabel('Time in Mega Years')
    #plt.ylabel('radius at plane intersect in au')
    plt.legend()
    #plt.title('minus')
    #plt.suptitle('Radius of Orbits on the Line of Intersection')
    #plt.suptitle('(a1=%s,e1=%s), (a2=%s,e2=%s)' % (a1 / au, e1, a2 / au, e2))

    #Difference Graph
    plt.figure()
    plt.subplot(121)

    plt.semilogy(T, abs(R1plus-R2plus), label='dRplus')
    plt.ylabel('diffenerence in radius intersect in m')
    plt.xlabel('Time in Mega Years')
    plt.legend()
    plt.title('plus')

    plt.subplot(122)
    plt.semilogy(T, abs(R1minus-R2minus) ,label='dRminus')
    # plt.xlabel('Time in Mega Years')
    # plt.ylabel('radius at plane intersect in au')
    plt.legend()
    plt.title('minus')
    plt.suptitle('(a1=%s,e1=%s), (a2=%s,e2=%s)' % (a1 / au, e1, a2 / au, e2))


    plt.show()
    return

def intersectfinder():
    Tmax=int(10*(10**6.))  #find all points in 100Myrs

    s1i=0
    s2i=0
    T=np.arange(0,Tmax,1e-1) #1e-1 year intervals
    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)

    dRplus=npr(a1,e1,(pi/2)+s1i+(wp1*(10**-6.)*T))-npr(a2, e2, (pi / 2) + s2i + (wp2*(10**-6.) * T))
    dRminus = npr(a1, e1, -(pi / 2) + s1i + (wp1*(10**-6.)*T))- npr(a2, e2, -(pi / 2) + (wp2*(10**-6.) * T))
    CrossDataplus=np.array([[0],[0],[0]])#,[0],[0]])
    CrossDataminus = np.array([[0], [0], [0]])
    for t in range(0,Tmax-1):
        if dRplus[t]*dRplus[t+1]<0: # if sign change
            CrossDataplus=np.concatenate((CrossDataplus,[[t],[((pi/2)+s1i+(wp1*(10**-6.)*t))%(2*pi)],[((pi / 2) + s2i + (wp2*(10**-6.) * t))%(2*pi)]]) ,axis=1)
        if dRminus[t]*dRminus[t+1]<0: # if sign change
            CrossDataminus=np.concatenate((CrossDataminus,[[t],[((-pi/2)+s1i+(wp1*(10**-6.)*t))%(2*pi)],[((-pi / 2) + s2i + (wp2*(10**-6.) * t))%(2*pi)]]) ,axis=1)
    CrossDataplus = np.delete(CrossDataplus, 0,1)
    CrossDataminus = np.delete(CrossDataminus, 0, 1)
    #print(CrossDataplus)
    return (CrossDataplus,CrossDataminus)

def saveintersects():
    (CrossDataplus, CrossDataminus)=intersectfinder()
    np.save('CrossDataplus.npy', CrossDataplus)#,CrossDataminus)
    np.save('CrossDataminus.npy', CrossDataminus)
    return

'''
saveintersects()

CrossDataplus=np.load('CrossDataplus.npy')
print(CrossDataplus)
'''





#intersectgraphs()
(CrossDataplus,CrossDataminus)=intersectfinder()
print(CrossDataplus)