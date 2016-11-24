from Old.OldOrbitCollisionGraph import *
i=1
for i in range (1,5):
    k=1
    i = i + 1
    for k in range (1,10):
        print('e1'),print(i/10),print('e2'),print(k/10)
        EllipseCrossImpact(2.2*au,i/10,0,2*au,k/10,2,2e30,2e10,2e10)
        plt.show(1)

        k=k+1