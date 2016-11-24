def EllipsePlot(a,e,s1):
    F = np.linspace(0, 2*pi, 100);
    R=npr(a,e,F-s1)
    plt.polar(F,R)
    return
def OldEllipseCrossImpact(a1,e1,s1,a2,e2,s2,Mstar,m1,m2):
    plt.figure(1)
    EllipsePlot(a1,e1,s1)
    EllipsePlot(a2,e2,s2)
    CollisionData=CollisionPoints(a1,e1,s1,a2,e2,s2)
    plt.scatter(CollisionData[0,:], CollisionData[1,:], c='r')
    #plt.scatter([0,pi/4,pi/2,pi],[0,10,15,20] ,c='b')
    [a3, e3, s3]=CollisionVars(a1, e1, s1, a2, e2, s2,Mstar,m1, m2)
    EllipsePlot(a3, e3, s3)
    return

def EllipseCross(a1,e1,s1,a2,e2,s2):
    plt.figure(1)
    EllipsePlot(a1,e1,s1)
    EllipsePlot(a2,e2,s2)
    CollisionData=CollisionPoints(a1,e1,s1,a2,e2,s2)
    plt.scatter(CollisionData[0,:], CollisionData[1,:], c='r')
    #plt.annotate(['Collision1','Collision2'],xy=(CollisionData[0,:], CollisionData[1,:]))