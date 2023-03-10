from pathlib import Path
import time
import csv
import numpy as np
import pandas as pd
import random
from collections import defaultdict 
from copy import deepcopy

def reduce(array):
    a = deepcopy(array)
    b = {a[0][i]: i+1 for i in range((a.shape[1]))} # map the first row to 1, 2, 3, ...
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            a[i][j] = b[a[i][j]]
    a = a[a[:, 0].argsort()]
    return (a)

def balance(n, k):
    # initial imbalance: return k*(n+1)/3
    # current imbalance: return the average of recriprocals of paired distances
    t = 0
    for i in range(1, n):
        for j in range(i+1, n+1):
            t += 1/abs(i-j)
    #print (n, k, t, k*t/n/(n-1)*2)
    return k*t/n/(n-1)*2

def rectangleImbalance(a):
    b = balance(a.shape[1], a.shape[0]) 
    imbalance = 0
    for u in range(1, a.shape[1]): 
        for v in range(u+1, a.shape[1]+1):
            local_imbalance = 0
            for row in range(a.shape[0]): 
                i, = np.where(a[row] == u)
                j, = np.where(a[row] == v)
                local_imbalance += 1/(abs(i-j))
                # compare the distances across every row, or something
            imbalance += abs(local_imbalance-b)
            #if abs(local_imbalance-b) > imbalance:
            #    imbalance = abs(local_imbalance-b)
    return imbalance

def orderImbalance(a):
    b = balance(a.shape[1], a.shape[0])
    imbalance = 0
    for u in range(1, a.shape[1]): 
        for v in range(u+1, a.shape[1]+1):
            local_imbalance = 0
            for row in range(a.shape[0]): 
                i, = np.where(a[row] == u)
                j, = np.where(a[row] == v)
                if i > j:
                    local_imbalance += 1
                else:
                    local_imbalance -= 1
            if a.shape[0] % 2 == 0:
                if local_imbalance != 0:
                    imbalance += 1
            else:
                if abs(local_imbalance) > 1:
                    imbalance += 1
                # compare the distances across every row, or something
    return imbalance

# Naive random Latin Square generation approaches

# Approach 0: Generate the latin square isomorphic to Zn.
def generateBoringLatinSquare(n):
    square = np.zeros([n, n], dtype=np.int8)
    for i in range(n):
        for j in range(n):         
            square[i][j] = ((i+j)%n)+1
    return square

# Approach 1: Generating random rows and appending them to the square if there is no conflict
# There is a slight optimization at the end, as the final row can be uniquely determined from the initial rows.
def generateRandomLatinSquare(n):
    #searches = 0  This is the number of randomly generated permutations that are made (which is effectively equivalent to the number of iterations)
    square = np.zeros([n, n], dtype=np.int8)
    columns = defaultdict(lambda: set())
    for i in range(n-1):
        #print (i)
        while True:
            isValid = 1
            a = np.random.permutation(n)+1
            for j in range(n):
                if a[j] in columns[j]: # duplicate found
                    break
            else:
                break
        # add the row to the square, and append individual values to the columns
        for j in range(n):
            columns[j].add(a[j])    
        square[i] = a
    l = set(range(1, n+1))
    a = []
    for j in range(n):
        a.append(list(l-columns[j])[0])
    square[-1] = a
    return square

# The average number of permutations generated for each value of n (from a sample size of 100) is as follows:
# {2: 1.0, 3: 3.86, 4: 13.71, 5: 60.06, 6: 347.23, 7: 2278.24, 8: 18711.61, 9: 148009.47}
# Obviously, this approach does not scale well to large numbers of latin squares. 

# Approach 2: Generating all latin squares of a fixed size and choosing one randomly. 
# This will not be done here, as even the n=7 case has 812851200 possibilities, which is only slightly more than the ~2300 iterations of the first approach.
# (According to OEIS, n=12 has merely 776966836171770144107444346734230682311065600000 possibilities.)

# Approach 3: The Jacobson+Matthew's approach to random generation
# Latin squares are represented by a nxnxn incidence matrix. Also, we would need a way to convert back:

def toIncidence(square):
    dim = square.shape
    if len(dim) != 2 or dim[0] != dim[1]:
        return "Not square"
    n = dim[0]
    incidenceMatrix = np.zeros([n, n, n])
    for i in range(n):
        for j in range(n):
            incidenceMatrix[square[i,j]-1,i,j] = 1
    return incidenceMatrix

def toLatinSquare(incidence):
    dim = incidence.shape
    if len(dim) != 3 or dim[0] != dim[1] or dim[0] != dim[2]:
        return "Not a valid incidence matrix"
    # confirm that all lines contain 1 one
    n = dim[0]
    square = np.zeros([n, n], dtype=np.int8)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if incidence[k, i, j] == 1:
                    square[i, j] = k+1
    return square

def newMatrix(incidence, perfect):
    # Initially, a random 0 tile is replaced with a -1
    dim = incidence.shape
    if len(dim) != 3 or dim[0] != dim[1] or dim[0] != dim[2]:
        return "Not a valid incidence matrix"
    n = dim[0]
    # Find a random (x0, y0, z0) with a value of 0
    if perfect == 1:
        while True:
            x0 = random.randint(0, n-1)
            y0 = random.randint(0, n-1) # TODO: swap in n or k idk
            z0 = random.randint(0, n-1)
            if incidence[x0, y0, z0] == 0:
                break
    else:
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if incidence[x][y][z] == -1:
                        x0, y0, z0 = x, y, z
    # Let x1 be the value of x such that (x1, y0, z0) = 1, and likewise for y1, z1
    xpos = []
    ypos = []
    zpos = []
    for x in range(0, n):
        if incidence[x, y0, z0] == 1:
            xpos.append(x)
    for y in range(0, n):
        if incidence[x0, y, z0] == 1:
            ypos.append(y)
    for z in range(0, n):
        if incidence[x0, y0, z] == 1:
            zpos.append(z)
    x1 = random.choice(xpos)
    y1 = random.choice(ypos)
    z1 = random.choice(zpos)
    incidence[x0, y0, z0] += 1
    incidence[x1, y1, z1] -= 1
    for x in [x0, x1]:
        for y in [y0, y1]:
            for z in [z0, z1]:
                if ((x==x0) and (y==y0) and (z==z0)) or ((x==x1) and (y==y1) and (z==z1)):
                    continue
                else:
                    incidence[x, y, z] = 1-incidence[x, y, z]
    if incidence[x1, y1, z1] == -1:
        perfect = 0
    else:
        perfect = 1
    return incidence,perfect

def shuffleMatrix(a): #shuffles two columns of a matrix
    b = deepcopy(a)
    c1 = random.randint(0, n-1)
    c2 = random.randint(0, n-1)
    temp = deepcopy(b[:, c1])
    b[:, c1] = b[:, c2]
    b[:, c2] = temp
    return b


values = {}
table_array = []
tables = {}

# read in all files in /tables
for file in Path('inverse').glob('*.csv'):
    numbers = file.stem[4:].split("_")
    i = int(numbers[1])
    j = int(numbers[0])
    tables[(i, j)] = np.genfromtxt(file, delimiter=',')
for i in tables:
    values[i] = rectangleImbalance(tables[i])
    
table = np.zeros([15, 15])
for i in values:
    table[i[0]][i[1]] = values[i]
for i in table:
    print (i)
DF = pd.DataFrame(table)
DF.to_csv("table.csv", index=False, header=False)

MAX_N = 12
for k in range(2, MAX_N+1):
    for n in range(k,MAX_N+1):
        if (n, k) in values and values[(n, k)] == 0:
            continue
        print ("Values: "+str(n)+" "+str(k))
        a = generateBoringLatinSquare(n)
        perfect = 1
        squares = defaultdict(lambda: 0)
        c = 0
        #print (a)
        minImbalance = 9999
        reset = 40
        for rng in range(120):
            # change matrix 
            a = toIncidence(a)
            reset -= 1
            while True:
                a, perfect = newMatrix(a, perfect)
                if perfect == 1:
                    break
            a = toLatinSquare(a)
            for i in range(35):
                b = shuffleMatrix(a)
                if rectangleImbalance(b[:k]) < rectangleImbalance(a[:k]):
                    a = b
                if rectangleImbalance(b[:k]) == 0:
                    table_array.append(reduce(a[:k]))
            r = rectangleImbalance(a[:k])
            if r <= minImbalance:
                minImbalance = r
                if type(r) == int:
                    values[(n, k)] = r
                else: 
                    values[(n, k)] = r[0]
                tables[(n, k)] = a[:k]
                if r == 0:
                    break
            if reset == 0:
                a = generateBoringLatinSquare(n)
                reset = 40
        DF = pd.DataFrame(tables[(n, k)])
        DF.to_csv("inversesquare/file"+str(k)+"_"+str(n)+".csv", index=False, header=False)

#for i in tables:
    # send each of the tables to a file named file_n_k.csv
    #w = csv.writer(fin, quoting=csv.QUOTE_ALL)


print (tables)