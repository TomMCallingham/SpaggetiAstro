from Functions import *
from PR.TprAnalytic import *
density=2000
Qa=620
Qb=(5e-6)*density
a=0.3
b=1.5
(a1,e1)=orbitalvalues(3)
(a2,e2)=orbitalvalues(10)
a1=0.999*a1
I2=1* ((2 * pi) / 360)
I1=0


#d2=1
def DispThres(D):
    Qd=(Qa*(D**-a))+(Qb*(D**b))
    return Qd
def DispShattering(D):
    Qs = (Qa * (D ** -a))
    return Qs
def flrfunc(vcol,D1,D2):
    Q = 0.5 * (vcol ** 2.)*((D2/D1)**3.)
    Qd=DispThres(D1)
    if Q<Qd: # cratering
        #print('Cratering')
        flr=1-0.5*(Q/Qd)
    if Q>=Qd: #shattering
        #print('Shattering')
        flr=0.5*((Qd/Q)**1.24)


    return flr

plt.show()