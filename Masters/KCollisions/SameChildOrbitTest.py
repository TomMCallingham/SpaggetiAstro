from Funcs import *
import numpy as np
import matplotlib.pyplot as plt
from KCollisions.KNewOrbit import *
Mstar=1.2e30
def ChildTest(a1,e1,s1,a2,e2,s2):
    F=np.linspace(0,2*pi,1000)
    rd1=rdot(a1,e1,F-s1,Mstar)
    td1=thetadot(a1,e1,F-s1,Mstar)
    rd2 = rdot(a2, e2, F - s2, Mstar)
    td2 = thetadot(a2, e2, F - s2, Mstar)
    Point=['a','b']
    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)  # Loading the ellipse intersection points
    print(CollisionData)
    plt.figure()

    plt.subplot(121)
    plt.plot(F / pi, (rd1 + rd2) / 2, label='average')
    plt.subplot(122)
    plt.plot(F / pi, (td1 + td2) / 2, label='average')

    for i in [0,1]:
        [a3, e3, s3, m3] = NewOrbit(a1, e1, s1, a2, e2, s2, Mstar, 1e8, 1e8,  CollisionData[0, i], CollisionData[1, i])
        rd3 = rdot(a3, e3, F - s3, Mstar)
        td3 = thetadot(a3, e3, F - s3, Mstar)

        plt.subplot(121)
        plt.plot(F / pi, rd3, label='child Point %s'%Point[i])
        plt.subplot(122)
        plt.plot(F / pi, td3, label='child Point %s' % Point[i])



    plt.subplot(121)
    plt.title('rdot')
    plt.scatter(CollisionData[0, 0] / pi, [0], label='Collision Angle a')
    plt.scatter(CollisionData[0, 1] / pi, [0], label='Collision Angle b')
    plt.legend()
    plt.subplot(122)
    plt.title('tdot')
    plt.scatter(CollisionData[0, 0] / pi, [0], label='Collision Angle a')
    plt.scatter(CollisionData[0, 1] / pi, [0], label='Collision Angle b')
    plt.legend()
    plt.show()
    return





ChildTest(2*au,0.8,0,2.1*au,0.6,2)
CollisionGraph(2*au,0.8,0,2.1*au,0.6,2,1.2e30,2e10,2e10)