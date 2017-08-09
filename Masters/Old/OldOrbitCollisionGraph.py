from KCollisions.KNewOrbit import *

def EllipseCrossImpact(a1,e1,s1,a2,e2,s2,Mstar,m1,m2):
    F = np.linspace(0, 2 * pi, 100)

    plt.figure()

#Plotting1
    R = npr(a1, e1, F - s1)
    plt.polar(F, R,label="Orbit1")
  #plotting 2
    # Plotting1
    R = npr(a2, e2, F - s2)
    plt.polar(F, R, label="Orbit2")
    CollisionData=CollisionPoints(a1,e1,s1,a2,e2,s2)
    plt.scatter(CollisionData[0,:], CollisionData[1,:], c='r')
    #plt.scatter([0,pi/4,pi/2,pi],[0,10,15,20] ,c='b')
    [a3, e3, s3,m3]=NewOrbit(a1, e1, s1, a2, e2, s2,Mstar,m1, m2,CollisionData[0,0],CollisionData[1,0])
    if e3<1 :
        R = npr(a3, e3, F - s3)
        plt.polar(F, R, label="Orbit3")
    else:
        print('Error: Orbit 3 Unbound!')

    [a4, e4, s4,m4] = NewOrbit(a1, e1, s1, a2, e2, s2, Mstar, m1, m2, CollisionData[0,1], CollisionData[1, 1])
    if e4 < 1:
        R = npr(a4, e4, F - s4)
        plt.polar(F, R, label="Orbit4")
    else:
        print('Error: Orbit 4 Unbound!')
    plt.legend()


au= 149597871e3
#EllipseCross(2*au,0,0,2*au,0.8,2)
EllipseCrossImpact(2*au,0.2,0,2.1*au,0.1,2,2e30,2e10,2e10),plt.show(1)


