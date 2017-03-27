import numpy as np

def discretise(v1,v2):
    #shift to v1=0, v2= v2-v1 and /h to get just numbers
    #do the same for masses?
    vI=v2-v1
    N=100
    J=10  #number of mass bins in each orbit
    T=1000

    h = vI / N #stepsize



    Data=np.zeros(M,N,T+1)  #Mass bins, Spacial position,timestep

    Space = np.linspace(v1, v2, N)

    #Seting up parent orbits

    #Collision prob matrix- does it depend on velocity?
    #Evolcing
    for t in range (0,T): #Evolving in time
        Data[:,:,t+1]=Data[:,:,t]
        for nparent in (0,N): #parent orbit 1 or 2

            for n2 in range (0,N):   # for each velocity bin
                #calc vimpact here!
                #collision loops  on each fragment size

                for jparent in range(1,J+1): #parent mass bins skipping the zeroth dead bin
                    for j2 in range (1,J+1):  #mass bin
                        #find the new orbit bin
                        n3=int(nparent+((j2/(j2+jparent))*(n2-nparent)))
                        #create prob*timestep create fragment vector

                        #take away used fragments!







