import numpy as np
import matplotlib.pyplot as plt
from KCollisions.KNewOrbit import *
from Funcs import *
G =6.674e-11# 2.982e-27
au= 1.496e11
Mstar=1.2e30
def eSolver(rd1,td1,rd2,td2,R):
    tdk = thetakepler(R, Mstar)

    #Remove R dimensions from r
    rd1=rd1/R
    rd2=rd2/R
    Dr=rd2-rd1
    Dt=td2-td1

    l=(Dt**4.)+((Dr**2.)*(Dt**2.))
    lA=(4*(Dt**3.)*td1)+((2*Dr*Dt)*((rd1*Dt)+(td1*Dr)))
    lB=(6*(Dt**2.)*(td1**2.))+((Dr**2.)*(td1**2.))+((Dt**2.)*(rd1**2.))+(4*rd1*Dr*td1*Dt)-(2*(tdk**2.)*(Dt**2.))
    lC=(4*Dt*(td1**3.))+((2*rd1*td1)*((Dr*td1)+(Dt*rd1)))-(4*(tdk**2.)*td1*Dt)

    Roots=np.roots([4*l,3*lA,2*lB,lC]) #the roots
    #print('Roots',Roots)

    #find emin
    F = (l * (Roots ** 4.)) + (lA * (Roots ** 3.)) + (lB * (Roots ** 2.)) + (lC * Roots)
    alpha=Roots[np.argmin(F)]
    emin=np.sqrt(((min(F) + ((td1 ** 4.) + ((rd1 ** 2.) * (td1 ** 2.)) - (2 * (tdk ** 2.) * (td1 ** 2.)))) / (tdk ** 4.)) + 1)
    #print( 'emin', emin)

    #print('emin',emin)
    rdmin=((1-alpha)*rd1+alpha*rd2)*R
    tdmin=(1-alpha)*td1+alpha*td2


    return [rdmin,tdmin]

def LambdaeData(a1,e1,a2,e2):




    N = 100# Choosing resolution badd word ah well
    L = np.linspace(0, 2 * pi, N+2)  # NOTE CHECK L def!! + or -
    L=L[1:N+1]
    R = np.zeros((2, N))  # radial collision points
    C = np.zeros((2, N))  # theta collision points

    for i in range(0, N):
        [[C[0, i], C[1, i]], [R[0, i], R[1, i]]] = CollisionPoints(a1, e1, 0, a2, e2, L[i])

    Rd1 = rdot(a1, e1, C[:, :], Mstar)  # another (2,N) array. Note taken out masses, no differnece!
    Td1 = thetadot(a1, e1, C[:, :], Mstar)
    Rd2 = rdot(a2, e2, C[:, :] - L, Mstar)
    Td2 = thetadot(a2, e2, C[:, :] - L, Mstar)
    #Tdk = np.sqrt((G * Mstar) / (np.power((R[:, :]), 3)))

    MineData=np.zeros((2,3,N))
    for j in range(0,2):
        for i in range(0, N):

            [MineData[j,1,i],MineData[j,2,i]]=eSolver(Rd1[j,i], Td1[j,i], Rd2[j,i], Td2[j,i], R[j,i]) # rdmin, tdmin
            MineData[j,0,i]=newe(MineData[j,1,i],MineData[j,2,i],R[j,i],Mstar) # min eccentriciy



    return MineData #0-emin1-rdmin,2-tdmin,


def LambdaeGraph(a1,e1,a2,e2):
    MineData=LambdaeData(a1,e1,a2,e2)
    N=np.size(MineData[0,0,:])
    L = np.linspace(0, 2 * pi, N + 2)  # NOTE CHECK L def!! + or -
    L = L[1:N + 1]
    plt.plot(L / pi, MineData[0, 0, :], label='BestEccentricity of Collision a')
    plt.plot(L / pi, MineData[1, 0, :], label='BestEccentricity of Collision b')
    plt.legend()
    plt.xlabel('Lambda/pi')
    plt.ylabel('Best Eccentricity')
    plt.show()
    return


LambdaeGraph(2*au,0.99,2*au,0.993)




