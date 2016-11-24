import numpy as np
A=np.array([[1,2,3],[4,9,6],[7,3,2]])
B=np.zeros((1,3))
C=np.zeros((1,3))
i=0
for i in range(0,3):

    C[0,i]=np.min(A[i,:])
    #B[0, i] = np.where(A[i,:] == C[0,i])
    B[0,i] = A[i,:].argmin()
    i +=1

print('B'),print(B)
print('C'), print(C)
print(A[0,2])