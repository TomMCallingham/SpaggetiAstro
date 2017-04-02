def fdistcalcOmega(Omega):
    Eafter=Epacket+((((m1*m2)/(m1+m2))*(((v2-v1))**2.))/Omega)
    print('Eafter/Ebefore',Eafter/Ebefore)
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
    plt.title('m1/m2=%s, Omega=%s,Eafter/Ebefore=%.3f, alpha1=%.2f'%(m1/m2,Omega,Eafter/Ebefore,alpha1))

    return
def OmegaGraphs():

    Omegas=[2.1,5,6,10,20]
    #print(np.size(Omegas))
    for i in range(0,np.size(Omegas)):
        print(i)
        fdistcalcOmega(Omegas[i])
    plt.show()
    return