import numpy as np
def n(N):
    n = np.zeros((2, N+1))  # number in row, total number
    # n[0:3]=[20, 11, 24] BETTER METHOD
    n[0, 0] = 2
    n[0, 1] = 1
    n[0, 2] = 2
    n[1, 0] = 2
    n[1, 1] = 3
    n[1, 2] = 5
    i = 3
    for i in range(3, N+1):
        n[0, i] = n[0, i - 1] * n[1, i - 2]
        n[1, i] = n[1, i - 1] + n[0, i]
        i += 1
    return n
print(n(9))