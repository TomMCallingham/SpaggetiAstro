import numpy as np
Mstar = 1.2e30
au = 1.496e11
year = 3.154e+7
#Setup
x='a'
a1=2*au
e1=0.995
s1=0
a2=2.1*au
e2=0.993
s2=1
dfSize = 1  #1 meter
def TPrecess(a,e):
    Tp=0.15*((1-(e**2.))/(1-(0.999**2.)))*((a/au)**2.5)
    wp = (2 * np.pi) / Tp
    return (Tp,wp)

(Tp1,wp1)=TPrecess(a1,e1)
(Tp2,wp2)=TPrecess(a2,e2)

Dwp=abs(wp1-wp2)
DTp=(2*np.pi)/Dwp
print('Tp1',Tp1)
print('Tp2',Tp2)
print('DTp',DTp)
print('Dwp',Dwp)
