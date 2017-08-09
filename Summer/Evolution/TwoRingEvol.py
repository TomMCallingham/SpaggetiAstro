from Evolution.RcolTContact import *
from Functions import *
from numba import jit, int32,float64

rmin=1
rmax=100
TotalVol=1e12 #((4*pi)/3)*(Rast**3.)
#Values
a1 = 2.5 * au
e1 = 0.998
a2 = 2.4 * au
e2 = 0.997
I=1* ((2 * pi) / 360)
Rbeam=1000

@jit
def InitialDistCalc(J):
    IniDist=np.zeros((2,J+1)) #number,rad
    IniDist[1,1:J+1]=np.linspace(rmin,rmax,J) #rad
    for f in range(1, J + 1):
        IniDist[0,f] = (  IniDist[1,f]) ** (-3.5) #Number in size range
        IniDist[0,0] += IniDist[0,f]*((IniDist[1,f])**3.)*((4*pi)/3) #Total Volume
    #print('TotalVol unscaled',Dist[0,0,0])
    IniDist[0, :] = ((TotalVol / IniDist[0, 0])) * IniDist[0, :]
    return IniDist
@jit
def StaticEvolution(s1,s2,x,T,tstep,J):
    #Setup
    # Data Matrix Setup
    Tnumber=int(T/tstep)
    Dist = np.zeros((3, J + 1, Tnumber + 1))  # P1P2C,Frag,Time
    IniDist = InitialDistCalc(J)
    #print('IniDist',IniDist[0, :])
    Dist[0, :, 0] = IniDist[0, :]
    Dist[1, :, 0] = IniDist[0, :]
    #CollisionRate Setup
    rcolwoA=RcolwoATcontact(s1,s2,x)[0]*tstep
    CrossSecMat = np.zeros((J + 1, J + 1))
    for i in range(1, J + 1):
        for j in range(1, J + 1):
            CrossSecMat[i, j] = (IniDist[1,i] + IniDist[1,j]) ** 2.
    '''
    print(rcolwoA/tstep)
    #print(CrossSecMat)
    print(IniDist[:,1])
    print('Thalf=',tstep/(IniDist[0,1]*rcolwoA*4))
    
    FirstCol=np.zeros((J,J))


    for f1 in range(1, J + 1):
        for f2 in range(1, J + 1):
            FirstCol[f1-1,f2-1]=np.log10(rcolwoA*CrossSecMat[f1,f2]*Dist[0,f1,0]*Dist[1,f2,0])
    print(FirstCol)
    '''

    #Loop

    for t in range(0,Tnumber): #Time Loop
        Dist[:,1:J+1,t+1]=Dist[:,1:J+1,t]

        for f1 in range(1,J+1):
            for f2 in range(1,J+1): # note mirrored problem!
                colnumber=rcolwoA*CrossSecMat[f1,f2]*Dist[0,f1,t]*Dist[1,f2,t] #Calc the number
                #print(colnumber)
                #Remove from parents
                Dist[0,f1,t+1]-=colnumber
                Dist[1,f2,t+1]-=colnumber
                #add to Child
                Dist[2,f1,t+1]+=colnumber
                Dist[2,f2,t+1]+=colnumber

    for f in range(1,J + 1):
        Dist[:,0,1:]+=Dist[:,f,1:]*((4*pi)/3)*(IniDist[1,f]**3.)


    return Dist

Dist=StaticEvolution(0,1,'a',1000,1e-1,10)
'''
print('Initial')
print(Dist[:,:,0])
print('Final')
print(Dist[:,:,np.size(Dist[0,0,:])-1])
TotalIn=Dist[0,0,0]+Dist[1,0,0]#+Dist[2,0,0]
TotalOut=Dist[0,0,Tnumber]+Dist[1,0,Tnumber]+Dist[2,0,Tnumber]
print('Totalin',TotalIn)
print('TotalOut',TotalOut)
print('Mass Change', TotalOut-TotalIn)
'''










