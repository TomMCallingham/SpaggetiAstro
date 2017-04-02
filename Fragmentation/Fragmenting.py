import numpy as np
import matplotlib.pyplot as plt
density=2700 #emily suggestioj
Qa=620
Qb=(5e-6)*density
a=0.5
b=1.5
vcol=1e5
d2=1
def DispThres(D):
    Qd=(Qa*(D**-a))+(Qb*(D**b))
    return Qd

def DispGraph():
    N=100
    D=np.linspace(-20,12,N)
    Q=np.zeros((N))
    for i in range(0,N):
        D[i]=10**D[i]
    Qd=DispThres(D)
    Q=0.5*(vcol**2.)#0.5*((d2/D)**3.)*(vcol**2.)
    #print('Q',Q)
    plt.figure()
    plt.loglog(D,Qd,label=' Dispersion')
    #plt.loglog([D[0],D[N-1]],[Q,Q],label='Spec Energy')
    plt.xlabel('Diameter D')
    plt.ylabel('Dispesion Thresh')
    plt.legend()
    plt.show()
    return


m1=1
m2=1
v1=10
v2=30
v3 = ((m1 * v1) + (m2 * v2)) / (m1 + m2)
Epacket=((m1+m2)*(v3**2.))/2
Ebefore=((m1*(v1**2.))+(m2*(v2**2.)))/2

print('Epacket/Ebefore', Epacket / Ebefore)
#Omega=100
#alpha1=2
def fdist1(v,v1,alpha1,K1):
    fd1=K1*((v-v1)**alpha1)
    return fd1
def fdist2(v,v2,alpha2,K2):
    fd2=K2*((v2-v)**alpha2)
    return fd2
def fdistcalc(Efrac):

    if Efrac<=(Epacket/Ebefore):
        print('Error, Energy below the packet')
    Eafter=Efrac*Ebefore

    #print('Eafter',Eafter)
    #print('Epacket',Epacket)
    #print('Eafter-Epacket',(Eafter-Epacket))
    Omega= (((m1*m2)/(m1+m2))*(((v2-v1))**2.))/(Eafter-Epacket)
    #print('Eafter/Ebefore',Eafter/Ebefore)
    alpha1=(+np.sqrt(1+(4*Omega))-5)/2
    alpha2=alpha1
    K1=m1*((1+alpha1)/((v3-v1)**(alpha1+1)))
    K2 = m2 * ((1 + alpha2) / ((v2 - v3) ** (alpha2 + 1)))
    v1space = np.linspace(v1, v3, 1000)
    v2space = np.linspace(v3, v2, 1000)
    fd1=fdist1(v1space,v1,alpha1,K1)
    fd2 = fdist2(v2space, v2, alpha2, K2)
    plt.figure()
    plt.plot((v1space-v1)/(v2-v1), fd1, label='fd1')
    plt.plot((v2space-v1)/(v2-v1), fd2, label='fd2')
    #plt.scatter(v3,0,'v3')
    plt.xlabel('position on velocity line')
    plt.ylabel('distribution')
    plt.legend()
    plt.title('m1/m2=%s, Eafter/Ebefore=%s alpha1=%.2f'%(m1/m2,Efrac,alpha1))

    return




def EfracGraphs():

    Efracs=np.linspace(Epacket/Ebefore+0.0001,0.9999,15)#[0.9001,0.91,0.92,0.93,0.94,0.95]
    #print(np.size(Omegas))
    for i in range(0,np.size(Efracs)):
        fdistcalc(Efracs[i])
    plt.show()
    return








EfracGraphs()