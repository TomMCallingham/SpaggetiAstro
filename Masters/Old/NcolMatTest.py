import numpy as np
J=1
N=11
MassFraction = np.zeros((J + 1, J + 1))
CrossSecMat = np.zeros((J + 1, J + 1))
for i in range(1, J + 1):
    for j in range(1, J + 1):
        MassFraction[i, j] = (i ** 3.) / ((i ** 3.) + (j ** 3.))

# collision bin calc

ncolmatunrounded = np.zeros((N, N, J + 1, J + 1))
for n1 in range(0, N-1 ):  # fragment 1
    for n2 in range(n1 + 1, N):
        ncolmatunrounded[n1, n2, :, :] = (n2 + n1) / 2# np.round( (n2 + n1)/2)
#ncolmatunrounded[:,:,1,1]=ncolmatunrounded[:,:,1,1]+np.transpose(ncolmatunrounded[:,:,1,1])
#ncolmatunrounded = ncolmatunrounded.astype(int)


ncolmatSymm = np.zeros((N, N, J + 1, J + 1))
for n1 in range(0, N-1 ):  # fragment 1
    ncolmatSymm[n1, N - (n1 + 1), 1, 1]= (N - 1) / 2
    for n2 in range(n1 + 1, N-(n1+1)):
        ncolmatSymm[n1, n2, :, :] =  np.round((n2 + n1) / 2)
        ncolmatSymm[N - (n2 + 1), N - (n1 + 1), :, :]= N - (ncolmatSymm[n1, n2, :, :] + 1)
#ncolmatunrounded[:,:,1,1]=ncolmatunrounded[:,:,1,1]+np.transpose(ncolmatunrounded[:,:,1,1])
ncolmatSymm = ncolmatSymm.astype(int)


# collision bin calc
ncolmatforcedup = np.zeros((N, N, J + 1, J + 1))
for n1 in range(0, N - 1):  # fragment 1
    for n2 in range(n1 + 1, N - (n1 )):
        ncolmatforcedup[n1, n2, :, :] = n1 + np.around((MassFraction[1, 1] * (n2 - n1)) + 0.001)
        ncolmatforcedup[N - (n2 + 1), N - (n1 + 1), :, :] = N - (ncolmatforcedup[n1, n2, :, :] + 1)
        '''
    for n2 in range(N - (n1 + 1), N):  # Velocity Bin Loop
        ncolmatforcedup[n1, n2, :, :] = n1 + np.around((MassFraction[1, 1] * (n2 - n1)) - 0.001)
        '''



ncolmatforcedup  = ncolmatforcedup .astype(int)

ncolmatforceddown = np.zeros((N, N, J + 1, J + 1))
for n1 in range(0, N - 1):  # fragment 1
    for n2 in range(n1 + 1, N - (n1 )):
        ncolmatforceddown[n1, n2, :, :] = n1 + np.around((MassFraction[1, 1] * (n2 - n1)) - 0.001)
        ncolmatforceddown[N - (n2 + 1), N - (n1 + 1), :, :] = N - (ncolmatforceddown[n1, n2, :, :] + 1)


ncolmatforceddown  = ncolmatforceddown .astype(int)

ncolmatforcedalt = np.zeros((N, N, J + 1, J + 1))
for n1 in range(0, N - 1):  # fragment 1
    for n2 in range(n1 + 1, N):
        ncolmatforcedalt[n1, n2, :, :] = n1 + np.around((MassFraction[1, 1] * (n2 - n1))+((-1)**(np.floor(n2/2)+np.floor(n1/2)))*0.001)  # rond to ceil



ncolmatforcedalt  = ncolmatforcedalt .astype(int)



# collision bin calc
ncolmatround = np.zeros((N, N, J + 1, J + 1))
for n1 in range(0, N - 1):  # fragment 1
    for n2 in range(n1 + 1, N):  # Velocity Bin Loop
        ncolmatround[n1, n2, :, :] = n1 + np.around(MassFraction[1, 1] * (n2 - n1))  # rond to ceil
ncolmatrounded = ncolmatround.astype(int)




#print(ncolmatunrounded[:, :, 1, 1])
#print('noclmatround',ncolmatround[:, :, 1, 1])
#print('ncolsymm',ncolmatSymm[0,:,1,1])
print('ncolforcedup',ncolmatforcedup [:,:, 1, 1])
print('ncolforceddown',ncolmatforceddown [:, :, 1, 1])
#print('ncolforcedalt',ncolmatforcedalt [0, :, 1, 1])
#print('ncol force-alt',ncolmatforcedup[:,:,1,1]-ncolmatforcedalt[:,:,1,1])
#print('ncol un-alt',ncolmatunrounded[:,:,1,1]-ncolmatforcedalt[:,:,1,1])