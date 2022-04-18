#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 20:04:59 2022

@author: ricardobortothopker
"""

#climate control Steve Eppinger DSM example
from complexFuncs import *
import numpy as np
import pandas as pd
# from scipy import optimize
# from geneticalgorithm import geneticalgorithm as ga
from simanneal import Annealer
import random
import matplotlib.pyplot as plt
from add_up_to import abstractions
import operator as op
from functools import reduce
random.seed(555)
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
A = np.array([[5  ,0  ,1.5,0  ,0.5],
              [0  ,2  ,0.5,0  ,1.5],
              [2.5,1.5,1  ,.5 ,0],
              [0  ,  0,1.5,1  ,0],
              [1.5,.5 ,0  ,0  ,3]]) 
L13 = abstractions(len(CC))
# L5_2 = level_comb(5,2)
class DSM_complex_problem(Annealer):
    def __init__(self, state, DSM,hist={}):
        self.DSM = DSM
        self.hist=hist
        self.count=0
        self.L = abstractions(len(DSM))
        # self.levelA = levelA
        super(DSM_complex_problem, self).__init__(state)  # important!
    def move(self):
        # initial_energy = self.energy()
        ei = self.energy()
        a = random.randint(0, 100)
        if a==1:
            a = random.randint(0, len(self.L) - 1)
            self.state[1]  = self.L[a]
        else:
            idx = self.state[0]
            c = random.randint(0, len(idx) - 1)
            b = random.randint(0, len(idx) - 1)
            self.state[0][c], self.state[0][b] = self.state[0][b], self.state[0][c]
        self.count+=1
        
        
        en = self.energy()
        self.hist[self.count] = {}
        self.hist[self.count]['energy'] = en
        self.hist[self.count]['idx'] = self.state[0]
        self.hist[self.count]['level'] = self.state[1]
        
    def energy(self):
        idx = self.state[0]
        levelA = self.state[1]
        return perc_effort_compexity(DSM_idx_arrange(self.DSM,idx),levelA)
a = random.randint(0, len(L13) - 1)
optimzer = DSM_complex_problem([list(range(len(CC))),L13[a]],CC)
# a = random.randint(0, len(L5_2) - 1)
# optimzer = DSM_complex_problem([list(range(len(A))),L5_2[a]],A,2)
# optimzer.set_schedule(optimzer.auto(minutes=0.1))
   
optimzer.copy_strategy = "slice"
optimzer.Tmax=25000
optimzer.Tmin=2
optimzer.updates=100
optimzer.steps=100000
optimzer.save_state_on_exit=True
state, e = optimzer.anneal()
steves = perc_effort_compexity(CC,[[0,1,2,3,4],[5,6,7],[8,9,10,11,12]])#Steve's solution 5189.53
print(f"Steve's solution: {steves}")
print(f"SA solution: {perc_effort_compexity(DSM_idx_arrange(optimzer.DSM,state[0]),state[1])}")
print(DSM_idx_arrange(optimzer.DSM,state[0]))
hist = optimzer.hist
df = pd.DataFrame(hist).T
df['count'] = df['level'].apply(lambda x:sum( [ len(xx) for xx in x])-len(x))
df['n_clusters'] = df['level'].apply(lambda x:len(x))
df = df.sort_values(by='energy')
df.boxplot(column='energy',by='count')
df.boxplot(column='energy',by='n_clusters')
df.boxplot(column='energy',by=['n_clusters','count'])
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom  # or / in Python 2
for i in df['n_clusters'].unique():
    df[df['n_clusters']==i].boxplot(column='energy',by='count')
    