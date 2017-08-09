@jit
def oldPREvolve(a,e,T,Rfrag):
    Timestep=1*year
    #print(Timestep)
    TempT=int(10**3.)
    ShowT=int(T/TempT)
    Data=np.zeros((3,ShowT+1)) #a,e
    TempData=np.zeros((3,TempT+1))
    Data[0,:]=np.arange(0, ShowT+1)*Timestep*ShowT
    Data[1:,0]=[a,e]
    L=Luminosity((np.arange(0,T+1))+0.5e9)
    #plt.plot(Data[0, :], L, label='L')
    econstant=((5*Q)/(8*np.pi*(c**2.)))*(3/(Rfrag*Density))
    aconstant=(Q/(4*np.pi*(c**2.)))*(3/(Rfrag*Density))
    for tshow in range(0,ShowT):
        TempData[1:,0]=Data[1:,tshow]

        for t in range(0,TempT):
            l=L[t+(tshow*TempT)]
            TempData[1, t + 1] = TempData[1, t] - Timestep*(aconstant * L[t+(tshow*TempT)] * (2 + (3 * (TempData[2, t] ** 2.)))) / (
            TempData[1, t] * ((1 - (TempData[2, t] ** 2.)) ** 1.5)) #a
            TempData[2,t+1]=TempData[2,t]-Timestep*((econstant*L[t+(tshow*TempT)]*TempData[2,t])/((TempData[1,t]**2.)*(np.sqrt(1-(TempData[2,t]**2.)))))#e
            if TempData[1, t + 1]<(10**-2)*au:
                TempData[1, t + 1]=(10**-2)*au
                TempData[2,t+1]=0
                break
        Data[1:,tshow+1]=TempData[1:,TempT]
        if Data[1, tshow + 1] < (10 ** -2) * au:
            Data[1, tshow + 1] = (10 ** -2) * au
            Data[2, tshow + 1] = 0
            break
    return Data