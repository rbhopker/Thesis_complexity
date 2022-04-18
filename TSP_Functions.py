#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 21:52:16 2021

@author: ricardobortothopker
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import pi,inf,exp,sqrt,factorial
from random import shuffle
import operator
from complexFuncs import complexity,DSM_rearrange

from scipy.spatial import ConvexHull, convex_hull_plot_2d

def distance(pt1,pt2):
    dx = pt1[0]-pt2[0]
    dy = pt1[1]-pt2[1]
    
    return (dx**2+dy**2)**.5
def calculate_euclidean_distance_matrix(xy):
    matrix = np.zeros([len(xy),len(xy)])
    for i in range(len(xy)):
        for j in range(i,len(xy)):
            dist = distance(xy[i],xy[j])
            matrix[i,j] = dist
            matrix[j,i] = dist
    return matrix
def total_distance(routine,dist_matrix):
    #calculates the total distance of a routine in the form [2,3,1,0,...]
    #dist a squared matrix with the distance between points
    dist = 0
    for i in range(len(routine)):
        dist +=dist_matrix[routine[i],routine[(i+1)%len(routine)]]
    return dist
        
    
    

def animatePlot(ax,fig,best_hist,it_hist,accep_hist):
    ax.clear()
    x_best,y_best,sols_best=list(zip(*best_hist))
    ax.plot(x_best,y_best,c='b')
    x_it,y_it,sols_it=list(zip(*it_hist))
    ax.scatter(x_it,y_it,c='k',label='iter')
    x_ac,y_ac,sols_sc=list(zip(*accep_hist))
    ax.scatter(x_ac,y_ac,c='r',label='accepted')
    plt.draw()
def SA_TSP(f,fneigh,x0,params,neighParams,f_params):
    """Peforms simulated annealing to find (optimum) solution for the Travel Salesperson"""
    initial_temp = f_params[0]
    final_temp = f_params[1]
    alpha = f_params[2]
    best_hist =[]
    # temp =[]
    current_temp = initial_temp
    neq = f_params[4] #number of rearrangements accepted at a given T
    nmax=neq*round(sqrt(len(x0)))
    nfrozen=f_params[3]  # if neq hasn't been reached at nfrozen successive
    schedule = f_params[5]
    live_plots=False
    counter = 0
    frozen=0;  # exit flag for SA
    naccept=1; # number of accepted configurations since last temperature change
    Tlast=1;   # counter index of last temperature change
    it_hist=[]
    accep_hist=[]
    # Start by initializing the current state with the initial state
    xi = x0
    sol = xi
    it_hist.append([initial_temp,f(xi,*params),sol])
    if live_plots:
        fig,ax=plt.subplots()
        # plt.ion()
        # fig.show()
    while current_temp >= final_temp and frozen<nfrozen:
        neighbor = fneigh(xi,*neighParams) 
        # tf=False
        # while not tf:
        #     neighbor = fneigh(xi)
        #     tf = constraints(neighbor,*params)
        #     tf=tf[0]
        
        it_hist.append([current_temp,f(neighbor,*params),sol])
        # ET.append(f(neighbor,*params))
        
        # Check if neighbor is best so far
        dE = f(neighbor,*params)-f(xi,*params)
        counter +=1
        if dE<0:
            xi=neighbor
            naccept+=1
            # frozen = 0
            accep_hist.append([current_temp,f(neighbor,*params),sol])

        # if the new solution is better, accept it

        # if the new solution is not better, accept it with a probability of e^(-cost/temp)
        else:
            if dE / current_temp>709:
                coef = 709
            else:
                coef =dE / current_temp
            if np.random.uniform(0, 1) < exp(-coef):
                xi = neighbor
                accep_hist.append([current_temp,f(neighbor,*params),sol])
        if f(xi,*params)<f(sol,*params):
            sol=xi
            frozen = 0
            # temp.append(current_temp)
            # hist.append(f(sol,*params))

        # decrement the temperature
        if (naccept<neq)&((counter-Tlast)<nmax):
            continue
        elif (naccept<neq)&((counter-Tlast)>=nmax):
            frozen=frozen+1
            Tlast=counter
            naccept=0
            # temp.append(current_temp)
            best_hist.append([current_temp,f(sol,*params),sol])
            if schedule ==1:
                current_temp -= alpha
            elif schedule ==2:
                if current_temp-0.00000001<current_temp*alpha:
                    current_temp = -1
                current_temp=current_temp*alpha
            if live_plots:
                animatePlot(ax,fig,best_hist,it_hist,accep_hist)
                # plt.draw()
        elif naccept==neq:
            Tlast=counter
            naccept=0
            # temp.append(current_temp)
            best_hist.append([current_temp,f(sol,*params),sol])
            if schedule ==1:
                current_temp -= alpha
            elif schedule ==2:
                if current_temp-0.00000001<current_temp*alpha:
                    current_temp =-1
                current_temp=current_temp*alpha
            if live_plots:
                animatePlot(ax,fig,best_hist,it_hist,accep_hist)
                # plt.draw()
            
            

    return sol,best_hist,it_hist,accep_hist
def plotResults(sol,xy,dist_matrix = None,string ='',save = False):
    if dist_matrix == None:
        dist_matrix = calculate_euclidean_distance_matrix(xy)
    fig,ax=plt.subplots()
    sortedXY = [xy[i] for i in sol]
    sortedXY.append(sortedXY[0])
    x_sorted,y_sorted=list(zip(*sortedXY))
    ax.plot(x_sorted,y_sorted)
    plt.xticks([])
    plt.yticks([])
    ax.scatter(x_sorted,y_sorted,c='r')
    for i in range(len(sol)):
        ax.annotate(str(sol[i]),(x_sorted[i],y_sorted[i]))
    ax.annotate('{}Total Distance: {:.2f}'.format(string,total_distance(sol,dist_matrix)),xy=(0.00, 1.02), xycoords='axes fraction')
    if save:
        plt.savefig('TSP Result {}.png'.format(string.split('\n')[0]), bbox_inches="tight")
    plt.show() 
def NeighborFunc(x,*options):
    # if 'percMutat' not in options:
    #     options['percMutat'] = .2
    lenT = len(x)-1
    nMut = round(options[0]*len(x))
    if nMut<2:
        nMut = 2
    to_mutate = np.random.choice(x[1:], nMut, replace=False, p=[1/lenT]*lenT)
    
    xx = np.array(x.copy())
    locX = []
    for i in to_mutate:
        locX.append(np.where(x==i)[0][0])
    shuffle(to_mutate)
    for i in range(len(to_mutate)):
        xx[locX[i]]=to_mutate[i]
    # while to_mutate!=np.array([]):
    #     first, to_mutate = to_mutate[0], to_mutate[1:]
    #     xx[x==v]=first
    
    return xx
def is_pareto(costs, where,return_mask = True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    if where =='TopLeft':
        where=1
    elif where =='TopRight':
        where =2
    elif where =='BottomLeft':
        where =3
    elif where =='BottomRight':
        where =4
    elif where not in [1,2,3,4]:
        raise NotImplementedError()
    #lt<
    #gt>
    if where == 3:
        # op1 = 1
        op2 = operator.lt
    elif where ==2:
        # op1 = 1
        op2 = operator.gt
    elif where ==4:
        x = costs.transpose()[0]*-1
        y = costs.transpose()[1]
        costs = np.array([x,y]).transpose()

        op2 = operator.lt
    elif where ==1:
        x = costs.transpose()[0]*-1
        y = costs.transpose()[1]
        costs = np.array([x,y]).transpose()
        op2 = operator.gt
    
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
    # while op1(next_point_index,len(costs)):
        # nondominated_point_mask = np.any(costs>costs[next_point_index], axis=1)
        nondominated_point_mask = np.any(op2(costs,costs[next_point_index]), axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient
def is_pareto_any(costs, return_mask = True):
    tf = np.zeros(len(costs), dtype=bool)
    for i in [1,2,3,4]:
        curTf = is_pareto(costs,i)
        tf = curTf | tf
    return tf
def dsm_tsp(sol,xyArr):
    leng = len(sol)
    dsm = np.zeros((leng,leng))
    tf = is_pareto_any(xyArr)
    tf = hull_vertices(xyArr)
    tf = np.invert(tf)
    for i in range(leng):
        s=0.5
        # if tf[i]:
        #     s=+.5
        #     dsm[sol[i],sol[i]]=1
        # if tf[(i+1)%leng]:
        #     s=+.5
        # dsm[sol[i],sol[(i+1)%leng]] =+s
        # dsm[sol[(i+1)%leng],sol[i]] =+s
        if tf[i]:
            s+=.5
            dsm[i,i]=2
        else: dsm[i,i]=1
            
        if tf[(i+1)%leng]:
            s+=.5
        dsm[i,(i+1)%leng] +=s
        dsm[(i+1)%leng,i] +=s
    # dsm = DSM_rearrange(dsm,sol)
    return dsm
def complexity_tsp_given_points_full(n_inner,n_outer,alpha_in = 2,alpha_out = 1, beta_in=2,beta_out=1):
    t_points = n_inner+n_outer
    C1 = n_inner*alpha_in+n_outer*alpha_out
    C2 = n_inner*beta_in+n_outer*beta_out
    leng = t_points
    dsm = np.zeros((leng,leng))
    for i in range(leng):
        dsm[i,i]=0
        dsm[i,(i+1)%leng] =1
        dsm[(i+1)%leng,i] =1
    SVD = np.linalg.svd(dsm, compute_uv=False) # SVD only in the architecture topology, no weight given to each interface
    C3 = sum(SVD)/leng
    return {'C1':C1,'C2':C2,'C3':C3,'C':C1+C2*C3}
def complexity_tsp_given_points(n_inner,n_outer,alpha_in = 2,alpha_out = 1, beta_in=2,beta_out=1):
    # t_points = n_inner+n_outer
    # C1 = n_inner*alpha_in+n_outer*alpha_out
    # C2 = n_inner*beta_in+n_outer*beta_out
    # leng = t_points
    # dsm = np.zeros((leng,leng))
    # for i in range(leng):
    #     dsm[i,i]=0
    #     dsm[i,(i+1)%leng] =1
    #     dsm[(i+1)%leng,i] =1
    # SVD = np.linalg.svd(dsm, compute_uv=False) # SVD only in the architecture topology, no weight given to each interface
    # C3 = sum(SVD)/leng
    C = complexity_tsp_given_points_full(n_inner,n_outer,alpha_in,alpha_out, beta_in,beta_out)
    return C['C']
def cal_tsp_complexity(xyArr):
    
    return complexity(dsm_tsp(list(range(len(xyArr))),xyArr))
def cal_tsp_complexityMult(xyArrTot):
    C = []
    for xyArr in xyArrTot:
        C.append(cal_tsp_complexity(xyArr))
    return C
def hull_vertices(xyArr):
    hull = ConvexHull(xyArr)
    tf = np.zeros(len(xyArr), dtype=bool)
    for ver in hull.vertices:
        tf[ver] = True
    return tf
def create_TSP_dumb(outer,total):
    ar = np.random.rand(2*total)
    ar = ar.reshape((total,2))
    if sum(hull_vertices(ar))==outer:
        return ar
    else:
        return create_TSP_dumb(outer,total)
def myround(x, prec=3, base=.05):
  return np.round(base * np.round(x/base),prec)    
def create_TSP_hull(outer,total,min_dist=0.05,tot_count=0):
    rng = np.random.default_rng()
    ar = np.random.rand(2*total)
    ar = ar.reshape((total,2))
    #ar = myround(ar,3,min_dist)
    cur_hull = hull_vertices(ar)
    def calc_dist(ar):
        return sum(calculate_euclidean_distance_matrix(ar)<min_dist)>1
    count = 0
    while sum(cur_hull)!=outer:
        if count ==100:
            tot_count+=1
            return create_TSP_hull(outer,total,min_dist*0.9,tot_count)
        if tot_count==5:
            return np.array([[0,0]])
        cur = sum(cur_hull)-outer
        # cur2 = calc_dist(ar)
        subAr = np.random.rand(2*abs(cur))
        subAr = subAr.reshape((abs(cur),2))
        if cur >0:
            tf = True
        else: tf = False
        tfArr0 = tf == cur_hull
        # tfArr1 = True == cur2
        # tfArr = tfArr0 | tfArr1
        loc = np.where(tfArr0)
        # subAr = np.random.rand(2*len(loc[0]))
        # subAr = subAr.reshape(len(loc[0]),2)
        # for i in range(abs(cur)):
        for i in range(abs(cur)):
            ar[loc[0][i]]=subAr[i]
        #ar = myround(ar,3,min_dist)
        rng.shuffle(ar)
        cur_hull = hull_vertices(ar)
        count +=1
    cur_dist = calc_dist(ar)
    count = 0
    while sum(cur_dist)>0:
        if count ==100:
            tot_count+=1
            return create_TSP_hull(outer,total,min_dist*0.9,tot_count)
        if tot_count==5:
            return np.array([[0,0]])
            
        cur = sum(cur_dist)
        tempAr = ar
        rng.shuffle(tempAr)
        subAr = np.random.rand(2*abs(cur))
        subAr = subAr.reshape((abs(cur),2))
        tfArr1 = True == cur_dist
        loc = np.where(tfArr1)
        for i in range(abs(cur)):
            tempAr[loc[0][i]]=subAr[i]
        if sum(hull_vertices(tempAr))==outer:
            ar = tempAr
            cur_dist = calc_dist(ar)
        count +=1
    return ar

def create_TSP_pareto(outer,total):
    rng = np.random.default_rng()
    ar = np.random.rand(2*total)
    ar = ar.reshape((total,2))
    pareto = is_pareto_any(ar)
    while sum(pareto)!=outer:
        cur = sum(pareto)-outer
        subAr = np.random.rand(2*abs(cur))
        subAr = subAr.reshape((abs(cur),2))
        if cur >0:
            tf = True
        else: tf = False
        loc = np.where(tf==pareto)
        for i in range(abs(cur)):
            ar[loc[0][i]]=subAr[i]
        rng.shuffle(ar)
        pareto = is_pareto_any(ar)     
    return ar
def create_TSP(outer,total,func):
    rng = np.random.default_rng()
    ar = np.random.rand(2*total)
    ar = ar.reshape((total,2))
    pareto = func(ar)
    while sum(pareto)!=outer:
        cur = sum(pareto)-outer
        subAr = np.random.rand(2*abs(cur))
        subAr = subAr.reshape((abs(cur),2))
        if cur >0:
            tf = True
        else: tf = False
        loc = np.where(tf==pareto)
        for i in range(abs(cur)):
            ar[loc[0][i]]=subAr[i]
        rng.shuffle(ar)
        pareto = func(ar)     
    return ar
def O_tsp(n):
    return factorial(n-1)/2
x=np.array([1])
y = pd.DataFrame()
y=1