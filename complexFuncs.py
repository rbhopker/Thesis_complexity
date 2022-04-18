# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 13:36:20 2021

@author: Ricardo Hopker
"""

import numpy as np
from itertools import permutations
from itertools import combinations
import pandas as pd
from functools import reduce
import matplotlib.pyplot as plt
# import pandas as pd
def complexity(DSM):
    DSM =np.array(DSM)
    DSM = np.nan_to_num(DSM)
    C1 = np.trace(DSM) #Component Complexity on the diagonal of the DSM
    np.fill_diagonal(DSM,0) # Clear Diagonal to calculate SVD and C2
    # C1 = len(DSM) 
    C2 = np.nansum(DSM) #Interface Complexity
    SVD = np.linalg.svd(DSM>0, compute_uv=False) # SVD only in the architecture topology, no weight given to each interface
    C3 = sum(SVD)/len(DSM) #Architecture (Topological) Complexity
    return {'C':C1+C2*C3,'C1':C1,'C2':C2,'C3':C3,'C2*C3':C2*C3,'SVD':SVD}

def DSM_rearrange(DSM,cols):
    outDSM = DSM
    rev = cols[::-1]
    outDSM[:,cols] = DSM[:,rev]
    outDSM[cols,:] = DSM[rev,:]
    # outDSM[cols,cols] = DSM[rev,rev]
    return outDSM
def DSM_idx_arrange(DSM,idx):
    outDSM = DSM
    outDSM = outDSM[:,idx]
    outDSM = outDSM[idx,:]
    return outDSM

def DSM_colapse_simple(DSM,idx):
    leng = len(DSM) - len(idx) + 1 # len of new DSM
    outDSM = np.zeros((leng,leng))
    NOTidx = list(range(len(DSM)))
    for j in idx:
        NOTidx.remove(j) # create an index of the matrix that will be unchanged
    MINidx = min(idx)
    MAXidx = max(idx)
    
    
    for i in NOTidx: # loop all lines of new DSM
        for j in NOTidx:
            if j>MINidx:
                colIDX = j-MAXidx+MINidx
            else: colIDX = j
            if i>MINidx:
                linIDX = i-MAXidx+MINidx
            else: linIDX = i
            # linIDX = i
            outDSM[linIDX,colIDX]=DSM[i,j]
    line = DSM[idx,:].sum(axis=0) # gets all the indexes to colapses and adds in line
    col = DSM[:,idx].sum(axis=1) # gets all the indexes to colapses and adds in column
    diag = complexity(DSM[idx][:,idx])['C'] # calculates the complexity of the colapse
    #newLine and newCol are the cross where the module colapses
    newLine = np.insert(line[NOTidx],MINidx,diag) # adds module's complexity to the diagonal in line 
    newCol = np.insert(col[NOTidx],MINidx,diag) # adds module's complexity to the diagonal in column 
    outDSM[min(idx),:] = newLine # assigns to the output the new line
    outDSM[:,min(idx)] = newCol # assigns to the output the new column
    return outDSM
def DSM_colapse_complex(DSM,idx): #complex colapsing. Maintains DSM complexity constant --- 0 based
    leng = len(DSM) - len(idx) + 1 # len of new DSM
    outDSM = np.zeros((leng,leng))
    NOTidx = list(range(len(DSM)))
    for j in idx:
        NOTidx.remove(j) # create an index of the matrix that will be unchanged
    MINidx = min(idx)
    MAXidx = max(idx)
    
    for i in NOTidx: # loop all lines of new DSM
        for j in NOTidx:
            if j>MINidx:
                colIDX = j-MAXidx+MINidx
            else: colIDX = j
            if i>MINidx:
                linIDX = i-MAXidx+MINidx
            else: linIDX = i
            # linIDX = i
            outDSM[linIDX,colIDX]=DSM[i,j]
    totalComplex = complexity(DSM) # system total complexity
    C2C3 = totalComplex['C2']*totalComplex['C3'] #Complexity due to architecture
    i = idx[0]
    outsideModuleComplex = complexity(DSM[NOTidx][:,NOTidx]) # calculates the complexity of the square matrix unchanged
    
    moduleComplex= complexity(DSM[idx][:,idx]) # calculates the complexity of the module
    remainingC2 = C2C3 - moduleComplex['C2']*moduleComplex['C3'] - outsideModuleComplex['C2']*outsideModuleComplex['C3']
    # line = DSM[idx,:].sum(axis=0) # gets all the indexes to colapses and adds in line
    line1 = DSM[NOTidx,:][:,idx]
    # col = DSM[:,idx].sum(axis=1) # gets all the indexes to colapses and adds in column
    col1 = DSM[:,NOTidx][idx,:]
    totSum = line1.sum()+col1.sum()
    # remainingC2 += totSum/(max(line1.shape)+max(col1.shape)) 
    # line1 = line1*remainingC2/totSum
    line1 = line1+line1*remainingC2/totSum
    line1=line1.sum(axis=1)
    # col1 = col1*remainingC2/totSum
    col1 = col1+col1*remainingC2/totSum
    col1=col1.sum(axis=0)
    diag = moduleComplex['C'] # gets the total complexity of the colapse
    #newLine and newCol are the cross where the module colapses
    line1 = np.insert(line1,MINidx,diag)
    col1 = np.insert(col1,MINidx,diag)
    
    outDSM[min(idx),:] = col1 # assigns to the output the new line
    outDSM[:,min(idx)] = line1 # assigns to the output the new column
    NewComplex = complexity(outDSM)
    
    #need to adjust outDSM comlexity by changing C2, since C1 and C3 cannot be changed
    toAdjustC2 = (NewComplex['C']- totalComplex['C'])/NewComplex['C3']
    
    line1[MINidx] = 0
    col1[MINidx] = 0
    totSum = line1.sum()+col1.sum()
    line1 = line1-line1*toAdjustC2/totSum
    col1 = col1-col1*toAdjustC2/totSum
    line1[MINidx] = diag
    col1[MINidx] = diag
    outDSM[min(idx),:] = col1 # assigns to the output the new line
    outDSM[:,min(idx)] = line1 # assigns to the output the new column

    
    return outDSM   
def DSM_colapse(DSM,idxs): #simple colapsing
    outDSM = DSM
    count = 0
    for idx in idxs:
        # diag = complexity(DSM[idx][:,idx])
        outDSM = DSM_colapse_simple(outDSM,np.array(idx)-count)
        # outDSM[count,count]=diag['C']
        count+=len(idx)-1
    return outDSM
def DSM_colapse2(DSM,idxs): #complex colapsing. Maintains DSM complexity constant
    outDSM = DSM
    count = 0
    # initComplex = CalculateComplexity(DSM)
    for idx in idxs:
        # if idx[0]==15:
        #     print(idx)
        # diag = CalculateComplexity(DSM[idx][:,idx])
        # print(np.array(idx)-(min(idx)-count))
        outDSM = DSM_colapse_complex(outDSM,np.array(idx)-count)
        # outDSM[count,count]=diag['C']
        # curComplex = CalculateComplexity(outDSM)
        count+=len(idx)-1
        # if abs(initComplex['C']-curComplex['C'])>0.0001:
        #     print(idx)
        #     print(initComplex['C']-curComplex['C'])
    return outDSM
def effort_complexity(complexity): #calculate effort based on given
    return 14.68*complexity**1.4775
def perc_effort_compexity(DSM,idx):#calculate effort to master perceived complexity 
    C = complexity(DSM)
    collapsed = DSM_colapse2(DSM,idx)
    count = 0
    effort = 0
    already_worked = 0
    for i in idx:
        if len(i)<=1:
            continue
        cur_lvl = i[0]-count
        count += len(i)-1
        
        effort += effort_complexity(collapsed[cur_lvl,cur_lvl])
        already_worked += collapsed[cur_lvl,cur_lvl]
    effort += effort_complexity(C['C']-already_worked)
    return effort
def level_comb(n,maxL):
    t = []
    t1 = []
    for i in range(2,n):
        cur = np.array(list(range(i)))
        for j in range(n-i+1):
            t.append(cur+j)
            t1.append([cur+j])
    combs={}
    combs[1] = t
    out =[]
    out.extend(t1)
    for i in range(2,maxL+1):
        LevelComb = list(combinations(t,i))
        valid = []
        for l in LevelComb:
            x=[]
            for xi in l:x.extend(xi)
            if len(x)==len(set(x)):
                valid.append(sorted(list(l),key=lambda x:x[0]))
        combs[i] = valid
        out.extend(valid)
    return out

# This code is contributed by mits

A = np.array([[5  ,0  ,1.5,0  ,0.5],
              [0  ,2  ,0.5,0  ,1.5],
              [2.5,1.5,1  ,.5 ,0],
              [0  ,  0,1.5,1  ,0],
              [1.5,.5 ,0  ,0  ,3]]) #pump controler exemple
# idxs = [[1,2],[3,4]]
# DSM_colapse2(A,idxs)
# perc_effort_compexity(A,idxs)
# possible_dsms = list(permutations(list(range(len(A)))))
# out = {}
# levels = level_comb(5,2)
# for perm in possible_dsms:
#     out[perm]={}
#     ident=0
#     for idx in levels:
#         out[perm][ident] = perc_effort_compexity(DSM_idx_arrange(A,perm),idx)
#         ident+=1
# df=pd.DataFrame(out)
# q = [1,1,1,1,2,2,2,3,3,2,2,3,2,3]
# fig,ax = plt.subplots()
# d1 =[]
# d2 =[]
# d3 =[]
# for col in df.columns:
#     d1.extend(df[col][0:4])
#     d2.extend(df[col][4:7])
#     d2.extend(df[col][9:11])
#     d2.extend([df[col][12]])
#     d3.extend(df[col][7:9])
#     d3.extend([df[col][11]])
#     d3.extend([df[col][13]])
#     ax.scatter(q,df[col],s=0.5,c='k')
# plt.boxplot([d1,d2,d3])

