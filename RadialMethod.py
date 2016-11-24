from Funcs import *
from NewOrbit import CollisionPoints

def RadialMethod(a1,e1,s1,a2,e2,s2,Mstar,m1,m2):
    [[ca, cb], [ra, rb]] =CollisionPoints(a1, e1, s1, a2, e2, s2)
    ra1dot=rdot(a1,e1,ca-s1,m1) #CHECK SIGN OF ANGLES
    ra2dot = rdot(a2, e2, ca - s2, m2)
    rb1dot = rdot(a1, e1, cb - s1, m1)
    rb2dot = rdot(a2, e2, cb - s2, m2)
    print('Orbit_a'), print(ra1dot),print(ra2dot)
    print('Orbit_b'), print(rb1dot),print(rb2dot)
    return

RadialMethod(2*au,0.993,0,2.1*au,0.99,1,1.2e30,2e10,2e10)