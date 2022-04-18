#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 00:25:09 2022

@author: ricardobortothopker
"""
import matplotlib.pyplot as plt
import numpy as np

class variable():
    def __init__(self,distribution,*args):
        self.dist = distribution
        self.args = args
    def evaluate(self):
        return self.dist(*self.args)
    def __call__(self):
        return self.evaluate()
    
class Monte_Carlo():
    def __call__(self):
        return self.evaluate()
    def simulation(self,n):
        self.results = [self.evaluate() for i in range(n)]
        return self.results
    def plot_pdf(self):
        fig = plt.figure()
        plt.hist(self.results,bins=25, density=True)
        plt.xlabel(self.x_label)
        plt.ylabel('PDF')
        plt.show()
    def plot_cdf(self):
        fig = plt.figure()
        plt.hist(self.results,bins=1000, density=True,cumulative=True, histtype='step',label = 'Cumulative Distribution')
        plt.xlabel(self.x_label)
        plt.ylabel('CDF')
        plt.vlines(np.array(self.results).mean(),0,1,color='r',label='Monte Carlo Average')
        plt.scatter(np.median(np.array(self.results)),0.5,c='r',label='Monte Carlo Median')
        plt.vlines(self.deterministic,0,1,color='k',linestyle='--',label='Deterministic')
        plt.legend(loc='lower right')
        plt.show()