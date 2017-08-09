import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from Funcs import *
from NewOrbit import *
def InclinationRot(I):
    RotMat=np.array([[cos(I),0,sin(I)],[0,1,0],[-sin(I),0,cos(I)]])



    return RotMat

def Orbitplot2D():
    a1 = 2
    e1 = 0.9
    s1 = 0
    a2 = 2.1
    e2 = 0.8
    s2 = 1

    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    x = 'a'
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]

    shift = (pi / 2) - C
    s1 = s1 + shift
    s2 = s2 + shift

    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure()
    ax = fig.add_subplot(111)
    N = 10000
    phi = np.linspace(0, 2 * np.pi, N)



    r = npr(a1, e1, phi - s1)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    plt.plot(x, y, label='orbit 1')

    r = npr(a2, e2, phi - s2)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    ax.plot(x, y, label='orbit 2')



    plt.scatter(CollisionData[0, 0], CollisionData[1, 0] / au, c='r')  # Highlighting the collision points
    plt.scatter(0, 0, c='b', marker='*')

    ax.set_xlabel('X')
    ax.set_xlim(-3, 3)
    ax.set_ylabel('Y')
    ax.set_ylim(-3, 3)
    ax.legend()

    plt.show()
    return

def Orbitplot3D():
    a1 = 2
    e1 = 0.9
    s1 = 0
    a2 = 2.1
    e2 = 0.8
    s2 = 1

    CollisionData = CollisionPoints(a1, e1, s1, a2, e2, s2)
    x='a'
    if x == 'a':
        R = CollisionData[1, 0]
        C = CollisionData[0, 0]
    elif x == 'b':
        R = CollisionData[1, 1]
        C = CollisionData[0, 1]
    #shift collision point to on the line
    shift=(-pi/2)-C
    s1 = s1 + shift
    s2 = s2 + shift

    #fig setup
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    N=10000
    phi = np.linspace(0, 2* np.pi, N)




    #parent 1
    r=npr(a1,e1,phi-s1)
    z = np.zeros((N))
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    ax.plot(x, y, z, label='Parent Ring 1',c='darkorange',zorder=0)

    #inclined parent 2
    r = npr(a2, e2, phi-s2)
    z = np.zeros((N))
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    RotMat=InclinationRot(-0.3)
    for i in range(0,N):
        X=np.array(([x[i],y[i],z[i]]))
        [x[i], y[i], z[i]]=RotMat@X

    ax.plot(x, y, z, label='Parent Ring 2',c='g',zorder=0)


    x = [0, 0]
    y = [-2, 3]
    z = [0, 0]
    ax.plot(x, y, z, label='Intersection Line', c='b',zorder=0)
    ax.scatter(0, 0, 0, marker='*', c='brown', zorder=9)
    ax.scatter(0, -CollisionData[1, 0], 0, c='r', zorder=10)

    ax.set_xlabel('X')
    ax.set_xlim(-1, 4)
    ax.set_ylabel('Y')
    ax.set_ylim(-2, 3)
    ax.set_zlabel('Z')
    ax.set_zlim(-1, 1.5)
    ax.legend()


    plt.show()
    return

#Orbitplot2D()

Orbitplot3D()