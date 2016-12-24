from NewOrbit import *
import numpy as np
au= 149597871e3
def RadialGraph(a1,e1,a2,e2):
    s1=0
    L = np.linspace(0, 2 * pi, 1000)
    R= np.zeros((2,1000))
    for i in range(1,1000):
        [[c1, c2], [R[0,i], R[1,i]]] = CollisionPoints(a1, e1, 0, a2, e2, -L[i])
    R=(1/au)*R
    plt.plot(L[1:1000],R[0,1:1000],label="Collision Point 1")
    plt.plot(L[1:1000], R[1, 1:1000], label="Collision Point 2")
    plt.plot([0,2*pi],[1,1])
    plt.legend() #attaching a legend to the graph
    plt.title('Radius of the Collision Points against Lambda')
    plt.xlabel('Lambda, Sepetation of the Orbits')
    plt.ylabel('Radius of the the Collision Points (au)')
    plt.show()
    return

RadialGraph(2*au,0.99,2.1*au,0.993)