import numpy as np
from Summer.Evolution.RcolTContact import *
from Summer.PR.TprAnalytic import *
from Summer.Fragmentation.FragFuncs import *
from numba import jit, int32,float64
I1=0
i2=10  #inlcination in degrees
I2=(i2/360)*(2*pi)
rsource1=3.1
rsource2=3
pmod=0.9
s1i=1
s2i=0
(a1,e1)=orbitalvalues(rsource1)
(a2,e2)=orbitalvaluesmod(rsource2,pmod)
rfrag = 1
rmincol=1e-6
Tmax=4e7
Tstep=1
TstepLarge=100

CrossDataFile='C:/Users/Tom/Documents/PycharmProjects/SpaggetiAstro/Summer/DataFiles/CrossData%s_%s_%s_%sMyr.npy' %(rsource1,rsource2,pmod,int(Tmax*(1e-6)))
MonoDistName='C:/Users/Tom/Documents/PycharmProjects/SpaggetiAstro/Summer/DataFiles/100RBMonoAccretion%s_%s_%s_%s_%sMyr.npy' %(rsource1,rsource2,pmod,i2,int(Tmax*(1e-6)))
ShrinkDistName='C:/Users/Tom/Documents/PycharmProjects/SpaggetiAstro/Summer/DataFiles/100RBShrinkMonoAccretion%s_%s_%s_%s_%sMyr.npy' %(rsource1,rsource2,pmod,i2,int(Tmax*(1e-6)))



def monoanalyticdelta(n1i,n2i,Rcol,t):
    delta=n2i-n1i
    n1=(n1i*delta*exp(-delta*Rcol*t))/((n1i*(1-exp(-delta*Rcol*t)))+delta)
    n2=delta+n1
    c=(n1i+n2i)-(n1+n2)
    dc=2*Rcol*n1*n2
    return(n1,n2,c,dc)

def monoanalytic(n1i,n2i,Rcol,t):
    n1=n1i/(1+(n1i*Rcol*t))
    n2=n1
    c=2*(n1i-n1)
    dc=2*Rcol*n1*n2
    return(n1,n2,c,dc)

def MonoAccretion(a1,e1,a2,e2):
    #Tstep=1
    Tnumber=int(Tmax/Tstep)

    Rast1=1e5
    Rast2=1e5
    nini1=(Rast1/rfrag)**3.
    nini2 = (Rast2 / rfrag) ** 3.

    print('Loading Data...')
    CrossData = np.load(CrossDataFile) #t, lambda=s1-s2,s1,s2,+/-1,a=0 b=1,R,C
    print('Creating MonoDist')
    AB=['a','b']
    Distributions=np.zeros((6,Tnumber)) #P1,P2,C,dC,Star,dStar
    M = np.size(CrossData[0, :])

    ContactData=np.zeros((4,M))# Rcol,Tcont,Tpr, flr3
    ChildData=np.zeros((4,M))# rd,td,a,e

    ParentVelocities=np.zeros((3,3,M))    #P1,P2,Mid
    Tendn = 0
    Distributions[0, 0] = nini1
    Distributions[1, 0] = nini2

    for m in range(0, M):
        print((m+1),'/',M)
        ContactData[0:2,m]=RcolwoATcontact(a1,e1,a2,e2,CrossData[2,m],CrossData[3,m],AB[int(CrossData[5,m])]) #Rcol, Tcontact
        ContactData[0,m]=ContactData[0,m]*4*(rfrag**2.)  #include cross sectional area!
        ParentVelocities[0, :, m] = np.array([rdot(a1, e1, CrossData[7,m]-CrossData[2,m]),CrossData[6,m]*thetadot(a1, e1, CrossData[7,m]-CrossData[2,m])*cos(I1),CrossData[6,m]*thetadot(a1, e1, CrossData[7,m]-CrossData[2,m])*sin(I1) ]) #p1 veloc
        ParentVelocities[1, :, m] = np.array(
            [rdot(a2, e2, CrossData[7,m]-CrossData[3, m]), CrossData[6, m] * thetadot(a2, e2, CrossData[7,m]-CrossData[3, m]) * cos(I2), #p2 veloc
             CrossData[6, m] * thetadot(a2, e2, CrossData[7,m]-CrossData[3, m]) * sin(I2)])
        ParentVelocities[2,:,m]=(ParentVelocities[0, :, m]+ParentVelocities[1, :, m])/2 #The Child Velocities
        ChildData[0:2,m]=np.array([ParentVelocities[2,0,m], np.sqrt((ParentVelocities[2,1,m] ** 2.) + (ParentVelocities[2,2,m] ** 2.))/CrossData[6,m]]) #Extracting rd, trd
        ChildData[2,m]=newa(ChildData[0,m],ChildData[1,m],CrossData[6,m])
        ChildData[3, m] = newe(ChildData[0,m],ChildData[1,m],CrossData[6,m])


        '''
        print('Lambda/pi',CrossData[1,m]/pi)
        print(AB[int(CrossData[5,m])])
        print('P1:rd,r*tdcosI,r*td*sinI', ParentVelocities[0, :, m])
        print('P2:rd,r*tdcosI,r*td*sinI', ParentVelocities[1, :, m])
        print('C:rd,r*tdcosI,r*td*sinI',ParentVelocities[2,:,m])
        print('collision rad',CrossData[6,m]/au)
        print('Child rd',ChildData[0,m])
        print('Child rtd', ChildData[1, m]*CrossData[6,m])
        print('new a ', ChildData[2,m]/au)
        print('new e ', ChildData[3, m])
        '''

        ContactData[2,m]=TPrAnalytic(ChildData[2,m],ChildData[3,m]) #Tpr
        ContactData[3,m]=((flrfunc(np.linalg.norm(ParentVelocities[0, :, m] - ParentVelocities[1, :, m]),2*rfrag,2*rfrag))**(1/3))*rfrag #largest surving
        if ContactData[3,m]<rmincol:
            print('Survivor smaller than 1 micron')
            ContactData[3, m] = rmincol


        #Collision Evolution Loop
        Tstartn=int((CrossData[0,m]-(ContactData[1,m]/2))/Tstep) #Find new start time
        Distributions[0,(Tendn+1):Tstartn+1]=Distributions[0,Tendn]  #Extend old distributions
        Distributions[1, (Tendn+1):Tstartn+1] = Distributions[1,Tendn]
        Distributions[2, (Tendn + 1):Tstartn + 1] += Distributions[2, Tendn] #Include Child
        Tendn = int((CrossData[0, m] + (ContactData[1, m] / 2)) / Tstep) #Find new end time
        Distributions[2, (Tstartn+1):(Tendn + 1)] += Distributions[2, Tstartn]  # Include Child
        #PR Drag eff
        topTPRn=int((ContactData[2,m]*ContactData[3,m])/Tstep) #find legnth of dist

        bottomTPRn=int((ContactData[2,m]*rmincol)/Tstep)
        PRAccrete=np.zeros((topTPRn-bottomTPRn))
        oldK=(rfrag**3.)/(2*(np.sqrt(ContactData[3,m]-rmincol)))
        '''
        print('Tpr', ContactData[2, m])
        print('Time Col', CrossData[0, m])
        print('Longest Tpr',topTPRn*Tstep)
        print('Fastest Tpr', bottomTPRn * Tstep)
        '''
        K=0
        for t in range(0,topTPRn-bottomTPRn):
            PRAccrete[t]=(t+bottomTPRn)**(-1/2)
            K+=PRAccrete[t]
        PRAccrete=   (1/K)*PRAccrete
        #PRAccrete=PRAccrete*K*(ContactData[2,m]**(-1/2))*(Tstep**2.)   #Check Tstep!

        
        #Collision Effects
        for t in range(Tstartn + 1, Tendn + 1):
            (n1, n2, c, dc) = monoanalytic(Distributions[0, Tstartn - 1], Distributions[1, Tstartn - 1],
                                           ContactData[0, m], (t-Tstartn) * Tstep)
            Distributions[0, t] = n1  # number in parent 1
            Distributions[1, t] = n2  # numer in parent 2
            Distributions[2, t] += c  # number in child need to fix!
            Distributions[3, t] = dc * Tstep  # ChildMassRate from collisoins

        # Effect of PR



        if (Tstartn+bottomTPRn+1)<Tnumber: #First Check PR will effect in time, remove weird errors



            if Tnumber > Tendn+topTPRn: #No Problem,  alwaysPR fits
                print('Fits')
                print('Range:',Tendn-Tstartn)
                for t in range(Tstartn + 1, Tendn + 1):

                    Distributions[2, (t + bottomTPRn):(t + topTPRn)] -= PRAccrete * Distributions[3, t] #Number in Child
                    Distributions[5, (t + bottomTPRn):(t + topTPRn)] += PRAccrete * Distributions[3, t] #Star Accretion Rate
            else: #PR accretion is longer than given time

                if Tnumber>Tstartn+topTPRn: #Some Fit some don't
                    #print('Some Fit')
                    #print('Range:',Tstartn + 1,Tnumber-topTPRn+1, 'Covering:', Tnumber-topTPRn - Tstartn)
                    for t in range(Tstartn + 1, Tnumber-topTPRn+1): #Fits
                        Distributions[2, (t + bottomTPRn):(t + topTPRn)] -= PRAccrete * Distributions[3, t]
                        Distributions[5, (t + bottomTPRn):(t + topTPRn)] += PRAccrete * Distributions[3, t]

                    #print('second half range',Tnumber-Tendn,'Tend',Tendn)
                    for t in range(Tnumber-Tendn, Tendn ): #Needs to be sliced
                        # Effect of PR
                        #print('Second half')
                        #print(t)
                        #print('start drop',(t + bottomTPRn))
                        Distributions[2, (t + bottomTPRn):] -= PRAccrete[0:(Tnumber-(t+bottomTPRn))] * Distributions[3, t]  # CHANGE
                        Distributions[5, (t + bottomTPRn):] += PRAccrete[0:(Tnumber-(t+bottomTPRn))] * Distributions[3, t]
                else: #None Fit
                    #print('None Fit')
                    #print('Range:', Tendn - Tstartn)
                    for t in range(Tstartn + 1, min(Tendn,Tnumber-bottomTPRn)):  # Needs to be sliced
                        # Effect of PR
                        Distributions[2, (t + bottomTPRn):] -= PRAccrete[:(Tnumber - (t + bottomTPRn))] * Distributions[
                            3, t]  # CHANGE
                        Distributions[5, (t + bottomTPRn):] += PRAccrete[:(Tnumber - (t + bottomTPRn))] * Distributions[
                            3, t]
        else: #If no PR material arrives
            print('no PR in time')
    Distributions[0, (Tendn + 1):] = Distributions[0, Tendn]  # Extend old distributions
    Distributions[1, (Tendn + 1):] = Distributions[1, Tendn]
    Distributions[2, (Tendn + 1):] += Distributions[2, Tendn]  # Include Child

    #Sum to Find Star Total
    Distributions[4, 0]=Distributions[5, 0] * Tstep
    for t in range(1, Tnumber):
        Distributions[4,t]=(Distributions[5,t]*Tstep)+ Distributions[4,t-1]

    Distributions=Distributions*((4*pi)/3)*density

    print('MonoDist Created')
    print('Largest Tpr', np.max(((ContactData[2,:].T)*ContactData[3,:])))

    return Distributions

def SaveMonoAccretionDist():
    MonoAccretionDist = MonoAccretion(a1,e1,a2,e2)
    print('saving')
    np.save(MonoDistName, MonoAccretionDist)
    print('saved')
    return
def ShrinkMonoAccretionDist():
    print('Loading Data...')
    MonoAccretionDist = np.load(MonoDistName)  #
    print('Creating Shrink Data...')

    largesmallratio=int(TstepLarge/Tstep)
    Tlargenumber=int(Tmax/TstepLarge)
    ShrinkDist=np.zeros((6,Tlargenumber))#P1,P2,C,dC,Star
    for tlarge in range(0,Tlargenumber):
        ShrinkDist[0:3,tlarge]=MonoAccretionDist[0:3,tlarge*largesmallratio] #P1,P2,Child
        ShrinkDist[4, tlarge] = MonoAccretionDist[4, tlarge * largesmallratio] #Star Total
        for rat in range(0,largesmallratio):
            ShrinkDist[3, tlarge]+=MonoAccretionDist[3,(tlarge*largesmallratio)+rat] #dchild
            ShrinkDist[5, tlarge] += MonoAccretionDist[5, (tlarge * largesmallratio) + rat] #dStar   #NEED TO DIVIDE?
    print('saving...')
    np.save(ShrinkDistName, ShrinkDist)
    print('saved')
    return

def MonoAccretionGraphs():
    print('Loading Data...')
    ShrinkDist = np.load(ShrinkDistName) ##P1,P2,C,dC,Star,dStar
    print('Plotting...')
    T=np.arange(0,Tmax,TstepLarge)
    TnumberLarge=int(Tmax/TstepLarge)
    TotalVol=np.zeros((TnumberLarge))
    for t in range(0,TnumberLarge):
        TotalVol[t]=ShrinkDist[0,t]+ShrinkDist[1,t]+ShrinkDist[2,t]#+ShrinkDist[4,t]

    plt.figure()
    plt.title('Mass in Parent Rings kg')
    plt.plot(T,ShrinkDist[0,:], label='P1')
    plt.plot(T, ShrinkDist[1, :], label='P2')
    plt.legend()
    plt.xlabel('Time, yrs')
    plt.figure()
    plt.title('Mass in Child kg')
    plt.plot(T, ShrinkDist[2, :], label='Child')
    plt.legend()
    plt.xlabel('Time, yrs')
    plt.figure()
    plt.title('Mass Accrete into Child From Collisions,g/s avaraged over 100yr periods')
    plt.semilogy(T, ShrinkDist[3, :]*(1000/(TstepLarge*year)), label='dChild')
    plt.legend()
    plt.xlabel('Time, yrs')
    plt.ylabel('grams a second')
    plt.figure()
    plt.title('PR driven Accretion onto Star, g/s avaraged over 100yr periods')
    plt.semilogy(T, ShrinkDist[5, :]*(1000/(TstepLarge*year)), label='Star Accretion Rate')
    plt.legend()
    plt.xlabel('Time, yrs')
    plt.ylabel('grams a second')
    plt.figure()
    plt.plot(T, ShrinkDist[4, :], label='TotalMass on Star')
    plt.ylabel('MAss')
    plt.legend()
    plt.xlabel('Time, yrs')
    plt.figure()
    plt.plot(T, TotalVol, label='Total Vol')
    plt.legend()
    plt.xlabel('Time, yrs')

    return
SaveMonoAccretionDist()
ShrinkMonoAccretionDist()
MonoAccretionGraphs()
plt.show()
