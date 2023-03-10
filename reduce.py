# reduces an ordered latin square via only tranposition and row shuffling (without column shuffling)
import numpy as np
a = np.array([[1, 6, 3, 5, 2, 4, 7],
    [2, 3, 4, 6, 7, 5, 1],
    [3, 7, 5, 4, 1, 6, 2],
    [4, 2, 6, 7, 3, 1, 5],
    [5, 1, 7, 3, 6, 2, 4],
    [6, 5, 1, 2, 4, 7, 3],
    [7, 4, 2, 1, 5, 3, 6]])
def reduce(a):
    b = {a[0][i]: i+1 for i in range(len(a[0]))} # map the first row to 1, 2, 3, ...
    for i in range(len(a)):
        for j in range(len(a[i])):
            a[i][j] = b[a[i][j]]
    a = a[a[:, 0].argsort()]
    print (a)
reduce(a)