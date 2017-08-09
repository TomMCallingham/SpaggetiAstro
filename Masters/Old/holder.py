[[ca, cb], [ra, rb]] = CollisionPoints(a1, e1, s1, a2, e2, s2)
Data1 = NewOrbit(a1, e1, s1, a2, e2, s2, Mstar, m1, m2, ca, ra)
print(Data1)
ra2dot = rdot(a2, e2, ca - s2, m2)
thetaa2dot = thetadot(a2, e2, ca - s2, m2)

m3 = m1 + m2
rdot3 = newmom(ra1dot, ra2dot, m1, m2, N1, N2)  # new rdot
thetadot3 = newmom(thetaa1dot, thetaa2dot, m1, m2, N1, N2)  # new thetadot
e3 = newe(rdot3, aData[N1, N2, 3], ra, aData[N1, N2, 1], Mstar)  # new e