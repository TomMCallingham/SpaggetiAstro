import numpy as np
from Functions import *
import matplotlib.pyplot as plt
Mstar = 1.2e30
au = 1.496e11
year = 3.154e+7
#Setup
rsource1=3
rsource2=10
pmod=0.9 #CHANGED
s1i=1
s2i=0
(a1,e1)=orbitalvalues(rsource1)
(a2,e2)=orbitalvaluesmod(rsource2,pmod)

Tmax=40*(10**6.)
tinterval = 1

CrossDataname='C:/Users/Tom/Documents/PycharmProjects/SpaggetiAstro/Summer/DataFiles/CrossData_rs%s_rs%s_pm%s_%sMyr.npy' %(rsource1,rsource2,pmod,int(Tmax*(1e-6)))
CrossID='(%s_%s_%s_%sMyr)'%(rsource1,rsource2,pmod,(Tmax*(1e-6)))
'''
#x='a'
a1=2.5*au
e1=0.998
#s1=0
a2=2.4*au #2.4
e2=0.997
#s2=1


p1=(1-e1)*a1
p2=(1-e2)*a2
q1=(1+e1)*a1
q2=(1+e2)*a2
if (p2-p1)*(q1-q2)>0:
    print('Rings Cross Twice Checked')
else:
    print('Rings dont cross twice')
'''
def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)*1e6 #in Yrs
    wp = (2 * np.pi) / Tp
    return (Tp,wp)

(Tp1,wp1)=TPrecess(a1,e1)
(Tp2,wp2)=TPrecess(a2,e2)





def intersectgraphs():
    T = np.arange(0, Tmax, 1000*tinterval)
    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)
    Tp1=Tp1*(10**-6) #convert to mega years
    Tp2 = Tp2 * (10 ** -6)

    R1plus=npr(a1,e1,(pi/2)-(s1i+wp1*T)) #changed to a minus
    R1minus = npr(a1, e1, -(pi / 2) -(s1i + wp1 * T))
    R2plus = npr(a2, e2, (pi / 2) - (s2i + wp2 * T))
    R2minus = npr(a2, e2, -(pi / 2) -( s2i + wp2 * T))

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
    plt.subplot(211)

    plt.semilogy(T, R1plus / au, label='R1plus')
    plt.semilogy(T, R2plus / au, label='R2plus')
    plt.ylabel('Radius at plane intersect in au')
    plt.title('Radius of Orbits on the Line of Intersection')
    #plt.xlabel('Time in Mega Years')
    plt.legend()
    plt.title('plus')

    plt.subplot(212)
    #plt.figure()
    plt.semilogy(T, R1minus / au, label='R1minus')
    plt.semilogy(T, R2minus / au, label='R2minus')
    plt.xlabel('Time in Mega Years')
    #plt.xlabel('Time in Mega Years')
    #plt.ylabel('radius at plane intersect in au')
    plt.legend()
    plt.title('minus')
    plt.suptitle('Radius of Orbits on the Line of Intersection')
    #plt.suptitle('(a1=%s,e1=%s), (a2=%s,e2=%s)' % (a1 / au, e1, a2 / au, e2))

    #Difference Graph
    '''
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

    '''
    #plt.show()
    return





def intersectfinder(): #note initiial has sli=1

      #find all points in Tmax years

    T=np.arange(0,Tmax,tinterval)
    (Tp1, wp1) = TPrecess(a1, e1)
    (Tp2, wp2) = TPrecess(a2, e2)

    dRplus=npr(a1,e1,(pi/2)-s1i-(wp1*T))-npr(a2, e2, (pi / 2) - s2i - (wp2 * T))   #changed from plue to minus!
    dRminus = npr(a1, e1, -(pi / 2) - s1i - (wp1*T))- npr(a2, e2, -(pi / 2)- s2i- (wp2 * T))
    CrossData=np.array([[0],[0],[0],[0],[0],[0],[0],[0]])  #t, lambda=s1-s2,s1,s2,+/-1,a=0 b=1,R,C
    print('searching for intersects...')
    for t in range(0,int((Tmax/tinterval)-1)):
        if dRplus[t]*dRplus[t+1]<0: # if sign change
            s1=(s1i+(wp1*t))%(2*pi) #+pi/2?
            s2= (s2i + (wp2 * t))%(2*pi)
            CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
            if abs(CollisionData[0,0]-(pi/2))<abs(CollisionData[0,1]-(pi/2)):  #its an a
                CrossData=np.concatenate((CrossData,[[t*tinterval],[(s1 - s2)% (2 * pi)],[s1],[s2],[+1],[0],[CollisionData[1,0]],[CollisionData[0,0]]]) ,axis=1)#[0]]) ,axis=1)
            else:#its a b
                CrossData = np.concatenate((CrossData, [[t * tinterval], [(s1 - s2)% (2 * pi)], [s1], [s2], [+1], [1],[CollisionData[1,1]],[CollisionData[0,1]]]), axis=1)  #[1]]), axis=1) #its a b
        if dRminus[t]*dRminus[t+1]<0: # if sign change
            s1 = ( s1i + (wp1  * t)) % (2 * pi)
            s2 = ( s2i + (wp2  * t)) % (2 * pi)
            CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
            if abs(CollisionData[0, 0] - ((3*pi) / 2)) < abs(CollisionData[0, 1] - ((3*pi) / 2)):  # its an a
                CrossData = np.concatenate((CrossData, [[t * tinterval], [(s1 - s2)% (2 * pi)], [s1], [s2], [-1], [0],[CollisionData[1,0]],[CollisionData[0,0]]]), axis=1)
            else:
                CrossData = np.concatenate((CrossData, [[t * tinterval], [(s1 - s2)% (2 * pi)], [s1], [s2], [-1], [1],[CollisionData[1,1]],[CollisionData[0,1]]]),
                                           axis=1)  #its a baxis=1)
    CrossData = np.delete(CrossData, 0,1)
    #print(CrossData)
    return (CrossData)



def saveintersects():
    (CrossData)=intersectfinder()
    print('saving')
    np.save(CrossDataname, CrossData)#,CrossDataminus)
    print('saved')
    return

def lambdapattern():
    CrossData = np.load(CrossDataname)

    print('plotting...')
    '''
    plt.figure()
    plt.scatter(CrossData[0,:],CrossData[1,:]/pi,c=((CrossData[4,:]+2)*(CrossData[5,:])))
    plt.xlabel('Time')
    plt.ylabel('lambda/pi in rad')
    plt.xlim([0, 20e6])
    plt.title('Value of Lambda at the Collision Points Coloured')

    plt.figure()
    plt.scatter(CrossData[0, :], CrossData[1, :] / pi, c=CrossData[4, :] )
    plt.xlabel('Time')
    plt.ylabel('lambda/pi in rad')
    plt.xlim([0, 20e6])
    plt.title('Value of Lambda at the Collision Points, plus minus')
    '''
    plt.figure()
    plt.scatter(CrossData[0, :], CrossData[1, :] / pi, c=CrossData[5, :])
    plt.xlabel('Time')
    plt.ylabel('lambda/pi in rad')
    plt.xlim([0, 20e6])
    plt.title('Lambda at the Collision Points (ab colour) %s'%CrossID)


    plt.figure()
    plt.scatter(CrossData[0,:],CrossData[7,:]/(pi),c=CrossData[4, :])
    plt.title('a or b check')
    '''
    
    plt.figure()
    plt.subplot(211)
    plt.scatter(CrossData[0, :], CrossData[2, :] / pi, c=CrossData[4, :])
    plt.xlabel('Time')
    plt.ylabel('s1/pi rad')
    plt.xlim([0,20e6])

    plt.subplot(212)
    plt.scatter(CrossData[0, :], CrossData[3, :] / pi, c=CrossData[4, :])
    plt.xlabel('Time in yrs')
    plt.ylabel('s2/pi in rad')
    plt.xlim([0, 20e6])

    plt.suptitle('Orientation of Rings at the Collision Points')
    #plt.show()
    '''
    plt.figure()
    plt.scatter(CrossData[0, :], np.log10(CrossData[6, :]/au))
    plt.xlabel('Time in yrs')
    plt.ylabel('R')
    plt.xlim([0, 20e6])


    return
#COMMANDS

#intersectgraphs()
#intersectfinder()
#saveintersects()
#CrossData=np.load(CrossDataname)
#print('no of contact:',np.size(CrossData[0,:]))
##print(np.transpose(CrossData))
#lambdapattern()
#print('Tp1/MYrs',Tp1), print('Tp2/Myrs',Tp2), print('wp1*Myrs',wp1), print('wp2*Myrs',wp2), print('wp1-wp2',wp1-wp2)
#Pp=(2*pi)/(wp1-wp2)
#print('PatternPeriod',Pp),print('no in 100',100/Pp),print('Pattern Density',4*(Pp/min(Tp1,Tp2))), print('expect in a Tmax:',(4/min(Tp1,Tp2))*Tmax)

#plt.show(), print('Plotted')