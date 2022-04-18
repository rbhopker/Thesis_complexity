# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 21:37:55 2021

@author: Ricardo Hopker
"""
import numpy as np
import pandas as pd
from complexFuncs import complexity, DSM_rearrange,DSM_idx_arrange,DSM_colapse_simple,DSM_colapse_complex,DSM_colapse,DSM_colapse2

url = r'C:\Users\Ricardo Hopker\OneDrive - Massachusetts Institute of Technology\Classes\Thesis'
matrix = pd.read_csv(url+r'\dsm_for_py.csv')
DSM = np.array(matrix)
matrix =np.array(matrix)
matrix = np.nan_to_num(matrix)
np.fill_diagonal(matrix,1) #simplest C1

C1 = np.trace(matrix)
np.fill_diagonal(matrix,0)
#simplest C2
C2 = np.nansum(matrix)
#C3
SVD = np.linalg.svd(matrix)
C3 = sum(SVD[1])/len(matrix)

print(C1+C2*C3)
modules =[0,4,8,15,19,21,25,32]
colapse=[list(range(0,4)),
         list(range(4,8)),
         list(range(8,15)),
         list(range(15,19)),
         list(range(19,21)),
         list(range(21,25)),
         list(range(25,32))]
# modules =[0,15,21,32]
s =[]
for i in range(len(modules)-1):
    s.append(sum(np.linalg.svd(matrix[modules[i]:modules[i+1]])[1]))
C3m =[]
C2m = []
C1m = []
Ctm = []
C=[]
for i in range(len(modules)-1):
    L = modules[i+1]-modules[i]
    A = matrix[modules[i]:modules[i+1],modules[i]:modules[i+1]]
    C3m.append(sum(np.linalg.svd(A)[1]/(L)))
    C2m.append(np.sum(A))
    C1m.append(L)
    Ctm.append(C1m[-1]+C2m[-1]*C3m[-1])
    C.append(complexity(A))
    
print(sum(Ctm))
      
DSM = np.array([[1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9],
    [1,2,3,4,5,6,7,8,9]])
idx = [3,4,5]
newDSM = DSM_colapse_simple(DSM,idx)

idxs = [[0,1,2],[3,4,5],[6,7,8]]
# newDSM2 = DSM_colapse(DSM,idxs)
test1 = DSM_colapse(matrix,colapse)
test1DF = pd.DataFrame(C)
tempDSM = np.array([[1,0,1],[0,2,2],[1,2,1]])
# DSM_colapse_two(tempDSM,[1,2])
# np.fill_diagonal(matrix,1)
comlexMatrix = complexity(matrix)
test11 = DSM_colapse_complex(matrix,list(range(0,4)))
comlexTest11 = complexity(test11)
test111 = DSM_colapse_simple(matrix,list(range(0,4)))
comlexTest111 = complexity(test111)
test2 = DSM_colapse2(matrix,colapse)
comlextest2 = complexity(test2)
test2DF = pd.DataFrame(C)      
A = np.array([[5  ,0  ,1.5,0  ,0.5],
              [0  ,2  ,0.5,0  ,1.5],
              [2.5,1.5,1  ,.5 ,0],
              [0  ,  0,1.5,1  ,0],
              [1.5,.5 ,0  ,0  ,3]]) #pump controler exemple
print(complexity(test1))
DSM = np.array([[1,2,3,4],[4,5,6,7],[7,8,9,0],[1,2,3,4]])
cols = [2,3]
newDSM = DSM_rearrange(DSM,cols)
DSM = np.array([[1,2,3,4,5],
                    [6,7,8,9,10],
                    [11,12,13,14,15],
                    [16,17,18,19,20],
                    [21,22,23,24,25]
                    ])
idx = [2,4,0,3,1]
print(DSM_idx_arrange(DSM,idx))