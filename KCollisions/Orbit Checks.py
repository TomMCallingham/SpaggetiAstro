import matplotlib.pyplot as plt
import numpy as np

from KCollisions.NewOrbit import *

au= 1#149597871e3
G = 2.982e-27#6.67408e-11

def OrbitChecks(a1, e1, a2, e2, Mstar):
    N = 1000  # Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N)  # NOTE CHECK L def!! + or -, basically s2, so opposite sign thatn i had before
    R = np.zeros((2, N))  # radial collision points
    C = np.zeros((2, N))  # theta collision points

    for i in range(1, N):
        [[C[0, i], C[1, i]], [R[0, i], R[1, i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    Rdot1 = rdot(a1, e1, C[:, :], Mstar)  # another (2,N) array. Note taken out masses, no differnece!
    Thetadot1 = thetadot(a1, e1, C[:, :], Mstar)
    Rdot2 = rdot(a2, e2, C[:, :] - L, Mstar)
    Thetadot2 = thetadot(a2, e2, C[:, :] - L, Mstar)

    # Eccentricity Check
    e1 = newe(Rdot1, Thetadot1, R, Mstar)
    e2 = newe(Rdot2, Thetadot2, R, Mstar)
    plt.figure()
    plt.subplot(121)
    plt.plot(L ,e1[0 ,:], label='e1 from a' )
    plt.plot(L, e2[0, :], label='e2 from a')
    plt.legend()
    plt.subplot(122)
    plt.plot(L ,e1[1 ,:], label='e1 from b' )
    plt.plot(L, e2[1, :], label='e2 from b')
    plt.legend()

    #EnergyCheck
    SE1=SEnergy(Rdot1, Thetadot1, R, Mstar)
    SE2 = SEnergy(Rdot2, Thetadot2, R, Mstar)
    plt.figure()
    plt.subplot(121)
    plt.plot(L, SE1[0, :], label='SE1 from a')
    plt.plot(L, SE2[0, :], label='SE2 from a')
    plt.legend()
    plt.subplot(122)
    plt.plot(L, SE1[1, :], label='SE1 from b')
    plt.plot(L, SE2[1, :], label='SE2 from b')
    plt.legend()

    plt.show()
    return


OrbitChecks(2*au,0.99,2.1*au,0.993,1.2e30)