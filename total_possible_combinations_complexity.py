#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 20:58:59 2022

@author: ricardobortothopker
"""


from simanneal import Annealer
import random
import matplotlib.pyplot as plt
from add_up_to import abstractions
import operator as op
from functools import reduce
import numpy as np
import random
from itertools import combinations, product, permutations
import time
from complexFuncs import perc_effort_compexity,DSM_idx_arrange,DSM_colapse2
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import string
from scipy.special import factorial

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom 

n = 13
lvls = abstractions(n)
def complex_comb(n,lvl):
    Num = n
    possibilities = 1
    denom = np.zeros(n)
    for cluster in lvl:
        r = len(cluster)
        possibilities *= ncr(Num,r)
        Num -= r
        denom[len(cluster)]+=1
        possibilities /= denom[len(cluster)]
    return possibilities
    
def total_complex_combinations(n):
    lvls = abstractions(n)
    total = 0
    for lvl in lvls:
        # possibilities = complex_comb(n,lvl)
        # Num = n
        # possibilities = 1
        # denom = np.zeros(n)
        # for cluster in lvl:
        #     r = len(cluster)
        #     possibilities *= ncr(Num,r)
        #     Num -= r
        #     denom[len(cluster)]+=1
        #     possibilities /= denom[len(cluster)]
        #have to divide possibilities by factorial of common cluster sizes
        total +=complex_comb(n,lvl)
        # total += possibilities
    return total
def possible_combinations(n,lvl):
    def combinations_from_list(lis,r):
        return list(combinations(lis,r))
    count = 0
    out = []    
    for cluster in lvl:
        r = len(cluster)
        cnr_possibilities = combinations_from_list(list(range(n-count)),r)
        out.append(cnr_possibilities)
        count+=r
    output = [[list(i)] for i in out[0]]
    for k in range(1,len(out)):
        output = [i+[list(j)] for i in output for j in out[k]]
    return output
def possible_combinations2(n,lvl):
    def combinations_from_list(lis,r):
        return list(combinations(lis,r))
    count = 0
    out = []   
    cluster = lvl[0]
    r = len(cluster)
    cnr_possibilities = combinations_from_list(list(range(n-count)),r)
    out = [set(frozenset(x)) for x in cnr_possibilities]
    for i in range(1,len(lvl)):
        cluster = lvl[i]
        r = len(cluster)
        for comb in out:
            remaining = set(range(n))
            for clus in list(comb):
                for line in clus:
                    remaining.remove(line)
            cnr_possibilities = combinations_from_list(list(remaining),r)
            
            
    # output = [[list(i)] for i in out[0]]
    # for k in range(1,len(out)):
    #     output = [i+[list(j)] for i in output for j in out[k]]
    return out
def all_combinations(n):
    lvls = abstractions(n)
    x = []
    for lvl in lvls:
        x = x+possible_combinations(n,lvl)
    return x
def all_combinations_from_lvls(n,lvls):
    x = dict()
    count = 0
    for lvl in lvls:
        combs = possible_combinations(n,lvl)
        x[count] = {'lvl':lvl,'combinations':combs,'flatten':flatten_comb_list(n,combs),
                    'clusters':len(lvl),'abstraction':len(flatten(lvl))-len(lvl),
                    'total':len(combs)}
        count+=1
    return x
def random_combination(n,lvl):
    def random_combination(iterable, r):
        "Random selection from itertools.combinations(iterable, r)"
        pool = tuple(iterable)
        n = len(pool)
        indices = sorted(random.sample(range(n), r))
        return tuple(pool[i] for i in indices)
    count = 0
    out =[]
    for cluster in lvl:
        r = len(cluster)
        cnr_possibilities = random_combination(list(range(n-count)),r)
        out.append(cnr_possibilities)
        count+=r
    return out
def flatten(list_of_lists):
    if len(list_of_lists) == 0:
        return list_of_lists
    if isinstance(list_of_lists[0], list):
        return flatten(list_of_lists[0]) + flatten(list_of_lists[1:])
    return list_of_lists[:1] + flatten(list_of_lists[1:])
# print(total_complex_combinations(n))
# all_combs = possible_combinations(n,lvls[10])
def flatten_comb(n,comb):
    x = list(range(n))
    out = []
    for i in comb:
        count = 0
        for j in i:
            out.append(x.pop(j-count))
            count += 1
    return out
def flatten_comb_list(n,combs):
    out = []
    for comb in combs:
        out.append(flatten_comb(n,comb))
    return out
def select_lvls(n,clusters,abstraction_level):
    lvls = abstractions(n)
    cluster = [x for x in lvls if len(x)==clusters]
    abstraction = [x for x in cluster if len(flatten(x))-clusters==abstraction_level]
    return abstraction
def show_matrix(DSM,idx,lvl,label):
    fig,ax = plt.subplots()
    ax.matshow(DSM_idx_arrange(DSM,idx))
    ticks = list(range(len(DSM)))
    label_sorted = label[idx]
    last = -.5
    for cluster in lvl:
        square = Rectangle((last,last), cluster[-1]+0.5-last, cluster[-1]+0.5-last,fc = None,ec="red",fill=False,lw=3)
        ax.add_patch(square)
        last = cluster[-1]+0.5
    ax.set_xticks(ticks,labels=label_sorted)
    ax.set_yticks(ticks,labels=label_sorted)
    fig.show()
    
CC = np.array([[1,1,0,1,1,0,0,0,0,0,0,0,0],
      [1,1,1,1,0,0,0,0,0,0,1,0,0],
      [0,1,1,1,0,0,0,0,0,0,0,1,1],
      [1,1,1,1,1,0,0,0,0,0,0,0,0],
      [1,0,0,1,1,1,0,0,1,0,0,0,0],
      [0,0,0,0,1,1,1,1,1,0,0,0,0],
      [0,0,0,0,0,1,1,1,0,0,0,0,0],
      [0,0,0,0,0,1,1,1,0,0,0,0,0],
      [0,0,0,0,1,1,0,0,1,1,1,0,0],
      [0,0,0,0,0,0,0,0,1,1,1,0,1],
      [0,1,0,0,0,0,0,0,1,1,1,1,1],
      [0,0,1,0,0,0,0,0,0,0,1,1,0],
      [0,0,1,0,0,0,0,0,0,1,1,0,1]])

CC = np.array([[1,1,1,1,1,0,0,0,0,0,0,0,0,], #unsorted
                        [1,1,1,0,0,0,0,0,0,0,1,1,0,],
                        [1,1,1,0,0,1,0,0,0,0,0,0,1,],
                        [1,0,0,1,1,0,0,0,0,0,0,0,0,],
                        [1,0,0,1,1,0,0,0,0,0,0,0,0,],
                        [0,0,1,0,0,1,1,0,0,0,0,0,1,],
                        [0,0,0,0,0,1,1,0,0,1,1,0,1,],
                        [0,0,0,0,0,0,0,1,0,1,1,0,0,],
                        [0,0,0,0,0,0,0,0,1,1,1,1,0,],
                        [0,0,0,0,0,0,1,1,1,1,0,0,1,],
                        [0,1,0,0,0,0,1,1,1,0,1,1,0,],
                        [0,1,0,0,0,0,0,0,1,0,1,1,0,],
                        [0,0,1,0,0,1,1,0,0,1,0,0,1,]])

cluster_3 = [x for x in lvls if len(x) in [2,3,4]]
abstraction_10 = [x for x in cluster_3 if len(flatten(x)) in [10,11,12,13]]
lvl = abstraction_10[6]
possible_combinations2(6,[[1,2],[1,2],[1,2]])
# meaningful_combs = all_combinations_from_lvls(n,abstraction_10)
# df = pd.DataFrame(meaningful_combs).T

# results = {}
# count = 0
# for i in range(len(df)):
#     lvl = df['lvl'].iloc[i]
#     combs = df['combinations'].iloc[i]
#     for comb in combs:
#         flat_comb = flatten_comb(n,comb)
#         current_DSM = DSM_idx_arrange(CC,flat_comb)
#         results[count]={}
#         results[count]['lvl'] = lvl
#         results[count]['Combination'] = comb
#         results[count]['Perceived Complexity'] = perc_effort_compexity(current_DSM,lvl)
#         results[count]['Flat'] = flat_comb
#         results[count]['clusters'] = len(comb)
#         results[count]['abstraction'] = len(flat_comb)-results[count]['clusters']
#         count+=1
# results = pd.DataFrame(results).T
# results.sort_values('Perceived Complexity',inplace=True)
# label = np.array(list(string.ascii_uppercase[:13]))
# show_matrix(CC,results['Flat'].iloc[0],results['lvl'].iloc[0],label)


# flatten_comb(n,meaningful_combs[0])
# result = []
# for comb in all_combs:
#     result.append(perc_effort_compexity(CC, comb))
# all_combs = all_combinations(8)
# for i in range(3,13):
#     start = time.time()
#     all_combs = all_combinations(i)
#     end = time.time()
#     print(f'n: %.0f; combinations: %.0f; time: %.3f'%(i,total_complex_combinations(i),end-start))


