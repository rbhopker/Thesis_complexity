#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 27 14:19:05 2022

@author: ricardobortothopker
"""

import pandas as pd
file_path = "/Users/ricardobortothopker/OneDrive - Massachusetts Institute of Technology/Classes/Thesis/excels/Points for TSP/"
file_name = "TSP summary_adjusted.xlsm"
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score,median_absolute_error,mean_absolute_percentage_error
from sklearn.metrics import explained_variance_score,max_error,mean_absolute_error,mean_squared_error,mean_squared_log_error
from TSP_Functions import complexity_tsp_given_points,O_tsp,complexity_tsp_given_points_full
from MonteCarlo import variable,Monte_Carlo
from random import triangular
import matplotlib.ticker as mtick
url = r'/Users/ricardobortothopker/OneDrive - Massachusetts Institute of Technology/Classes/Thesis/excels/Points for TSP/'
name = url+r'point count for TSP3.csv'

def get_ids():
    cols = {}
    for i in range(0,30,1):
        cols[i] = f"id {i}"
    cols[29] = "id 30"
    return cols
def load_results(file=file_path+file_name):
    df = pd.read_excel(file,sheet_name="TSP summary",skiprows=[0,1],nrows = 30).T
    rows = ['id','bestSol','sol','C','C1','C2','C3','n-nodes','n-points','sorted sol']
    count=1
    for i in range(len(rows),len(df),4):
        rows.extend([f'solution_{count}',f'check_{count}',f'cost_{count}',f'time_{count}'])
        count+=1
    df['row_names'] = rows
    df = df.set_index('row_names')
    cols = get_ids()
    df = df.rename(columns=cols)
    def get_id_values(df,ident):
        df_temp =df[ident].iloc[10:]
        rows_sol = [f'solution_{j}' for j in range(1,int(len(df_temp)/4+1))]
        rows_check = [f'check_{j}' for j in range(1,int(len(df_temp)/4+1))]
        rows_cost = [f'cost_{j}' for j in range(1,int(len(df_temp)/4+1))]
        rows_time = [f'time_{j}' for j in range(1,int(len(df_temp)/4+1))]
        df_sol = pd.DataFrame([df_temp.loc[rows_sol]]).T.reset_index().drop(columns=['row_names']).rename(columns={ident:'solution'})
        df_check = pd.DataFrame([df_temp.loc[rows_check]]).T.reset_index().drop(columns=['row_names']).rename(columns={ident:'check'})
        df_cost = pd.DataFrame([df_temp.loc[rows_cost]]).T.reset_index().drop(columns=['row_names']).rename(columns={ident:'cost'})
        df_time = pd.DataFrame([df_temp.loc[rows_time]]).T.reset_index().drop(columns=['row_names']).rename(columns={ident:'time'})
        df_out = pd.concat([df_sol,df_check,df_cost,df_time],axis=1)
        df_out['C'] = [df[ident].loc['C']]*len(df_out)
        df_out['id'] = [ident]*len(df_out)
        df_out['normalized_cost'] = df_out['cost']/df[ident].loc['bestSol']
        df_out = df_out.dropna(thresh=5).reset_index()
        # df_out['time/C'] = df_out['time']/df_out['C']
        return df_out
    all_df = {}
    unified_df = get_id_values(df,'id 0')
    for i in range(1,30,1):
        ident = f'id {i}'
        if i == 29:
            ident = 'id 30'
        all_df[ident] = get_id_values(df,ident)
        unified_df = pd.concat([unified_df,all_df[ident]])
    unified_df = unified_df.reset_index(drop=True)
    unified_df = unified_df.rename(columns={'index':'Test Id'})
    return unified_df
def regression_metric(y,ypred):
    out={}
    out['r2'] = r2_score(y, ypred)
    out['EVS'] = explained_variance_score(y, ypred)
    out['max_error'] = max_error(y, ypred)
    out['MAE'] = mean_absolute_error(y, ypred)
    out['MSE'] = mean_squared_error(y, ypred)
    out['MSLE'] = mean_squared_log_error(y, ypred)
    out['MedianAE'] = median_absolute_error(y, ypred)
    out['MAPE'] = mean_absolute_percentage_error(y, ypred)
    out['residual'] = y-ypred
    return out
def Custom_Fit(x,y,func):
    popt, pcov = curve_fit(func, x, y)
    # popt, pcov = curve_fit(func, x, y,epsfcn=)
    reg = {}
    reg['metrics'] = regression_metric(y,func(x, *popt))
    
    reg['coef']={}
    reg['coef']['a'] = popt[0]
    reg['coef']['b'] = popt[1]
    reg['coef']['d'] = popt[2]
    reg['cov'] = pcov
    return reg
def plot_and_fit(df, x,y,func,func_label = 'fit: a=%5.3f, b=%5.3f',xlabel='Complexity (C)',ylabel = 'Time (s)',ylim=None):
    fig,ax = plt.subplots()
    ax.scatter(x,y , c='b', label='data')
    popt, pcov = curve_fit(func, x, y)
    ax.plot(df['C'].unique(), func(df['C'].unique(), *popt), 'r-',
             label=func_label % tuple(popt))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend()
    if not ylim == None:
        ax.set_ylim(ylim)
    reg = {}
    reg['metrics'] = regression_metric(y,func(x, *popt))
    r2 = reg['metrics']['r2']
    ax.text(0,1.05,'r^2=%5.3f'%r2, horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes)
    fig.show()
    reg['coef']={}
    reg['coef']['a'] = popt[0]
    reg['coef']['b'] = popt[1]
    try:
        reg['coef']['d'] = popt[2]
    except:
        pass
    reg['cov'] = pcov
    return reg

def fit_func_1(x,a,b):
    return a*(x**b) 
def fit_func_2(x,a,b,c):
    return a*(x**b)+c
def analysis_with_plot(df):
    def plot_given_df(df):
        x = df['C']
        y = df['time']
        # results_1 = plot_and_fit(df,x,y,fit_func_2,func_label = 'fit: a=%5.3f, b=%5.3f, d=%5.3f',xlabel='Complexity (C)',ylabel = 'Time (s)',ylim=[0,250])

        results_1 = plot_and_fit(df,x,y,fit_func_1,func_label = 'fit: a=%5.3f, b=%5.3f',xlabel='Complexity (C)',ylabel = 'Time (s)',ylim=[0,200])
        y = df['normalized_cost']
        results_3 = plot_and_fit(df,x,y,fit_func_2,func_label = 'fit: a=%5.6f, b=%5.3f, d=%5.3f',xlabel='Complexity (C)',ylabel = 'n-optimium values',ylim=[1,1.2])
        y = df['tn']
        results_5 = plot_and_fit(df,x,y,fit_func_2,func_label = 'fit: a=%5.6f, b=%5.3f, d=%5.3f',xlabel='Complexity (C)',ylabel = 'n-optimium values * time',ylim=[0,250])

        results_1['tests'] = t_p_test(results_1,df)
        results_3['tests'] = t_p_test(results_3,df)
        results_5['tests'] = t_p_test(results_5,df)
        return [results_1,results_3,results_5]
    df['tn']=df["time"] * df["normalized_cost"]
    tf = np.invert(np.isnan(df['time']))
    df = df[tf]
    results_1,results_3,results_5 = plot_given_df(df)
    
    df2 = df[tf].groupby('C').mean()
    df2=df2.reset_index()
    results_2,results_4,results_6 = plot_given_df(df2)
    df_temp = df[['C','id']][tf].groupby('id').mean()
    df_temp['count'] = df[['C','id']][tf].groupby('id').count()
    df_temp['time'] = df[['id','time']][tf].groupby('id').mean()
    df_temp['timeSTD'] = df[['id','time']][tf].groupby('id').std()
    df_temp['cost'] = df[['id','cost']][tf].groupby('id').mean()
    df_temp['costSTD'] = df[['id','cost']][tf].groupby('id').std()
    df_temp['normalized_cost'] = df[['id','normalized_cost']][tf].groupby('id').mean()
    df_temp['normalized_costSTD'] = df[['id','normalized_cost']][tf].groupby('id').std()
    df_temp['tn'] = df[['id','tn']][tf].groupby('id').mean()
    df_temp['tnSTD'] = df[['id','tn']][tf].groupby('id').std()
    df3=df_temp.copy()
    for i in range(len(df_temp)):
        cur = df_temp.iloc[i]
        count = cur['count']
        for j in range(int(count)-1):
            df3= df3.append(cur)
    df3 = df3.sort_values(by='C')
    results_7,results_8,results_9 = plot_given_df(df3)
    
    
    results = pd.DataFrame({'Full -> t = aC^b':results_1,'Full -> n_opt = aC^b+d':results_3,'Full -> n_opt*t = aC^b+d':results_5,
                            'Av. -> t = aC^b':results_2,'Av. -> n_opt = aC^b+d':results_4,'Av. -> n_opt*t = aC^b+d':results_6,
                           'Weighted Av. -> t = aC^b':results_7,'Weighted Av. -> n_opt = aC^b+d':results_8,'Weighted Av. -> n_opt*t = aC^b+d':results_9})
    fig,ax = plt.subplots(3,3,sharex='col')
    c1 = 0
    c2 = 0
    for i in results.columns:
        ax[c1][c2].scatter(list(range(len(results[i]['metrics']['residual']))),results[i]['metrics']['residual'])
        ax[c1][c2].hlines(0,0,len(results[i]['metrics']['residual']),color='r')
        ax[c1][c2].set_title(i)
        c1 +=1
        if c1 >2:
            c2+=1
            c1=0
            
    fig.show()
    return results
def t_p_test(metric,df):
    std = np.sqrt(np.diag(metric['cov']))
    try:
        params = np.array([metric['coef']['a'],metric['coef']['b'],metric['coef']['d']])
    except:
        params = np.array([metric['coef']['a'],metric['coef']['b']])
    ts_b = params/ std
    p_values =[2*(1-stats.t.cdf(np.abs(i),(len(df)-1))) for i in ts_b]
    
    inter = [stats.t.interval(alpha=0.7, df=len(df)-1,loc = params[i],scale=std[i])
             for i in range(len(ts_b))]
    try:
        return pd.DataFrame({'value':[metric['coef']['a'],metric['coef']['b'],metric['coef']['d']],
                             'standard error':std,'t-test':ts_b,'p-value':p_values,'interval':inter},index=['a','b','d'])
    except:
        return pd.DataFrame({'value':[metric['coef']['a'],metric['coef']['b']],
                             'standard error':std,'t-test':ts_b,'p-value':p_values,'interval':inter},index=['a','b'])
def analysis(df):
    df['tn']=df["time"] * df["normalized_cost"]
    tf = np.invert(np.isnan(df['time']))
    df2 = df[tf].groupby('C').mean()
    x = df2.index.to_numpy()
    y = df2['tn'].to_numpy()
    metric = Custom_Fit(x, y, fit_func_2)
    metric['tests'] = t_p_test(metric,df2)
    return metric
def re_calculate_complexity(df,a_in,a_out,b_in,b_out):
    out = []
    df_tsp = pd.read_csv(name)
    df_out = df.copy()
    for i in range(len(df_tsp)):
        inner = df_tsp['interior points'].iloc[i]
        outer = df_tsp['outside points'].iloc[i]
        out.append(complexity_tsp_given_points(inner,outer,
                                               a_in,a_out,
                                               b_in,b_out))
    df_tsp['New_C'] = out
    df_out['C']= df_out['id'].map(df_tsp[['id','New_C']].set_index('id').squeeze())
    return df_out
def plot_contour_ratio(A,B,Z,over_1=True):
    fig, ax = plt.subplots()
    
    CS = ax.contour(A, B, Z,levels = np.arange(0.8,2,0.1))
    manual_locations = [
        (1, 1), (1.1, 1.1), (1.5, 1.5), (1.2, 1.2), (1.3, 1.2), (1.3, 1.3),(1.4,1.4)]
    ax.clabel(CS, inline=True, fontsize=10)
    ax.set_title('Exponent coefficient ''b', fontsize=20)
    ax.set_ylabel(r'$\frac{\beta_{inner}}{\beta_{outer}}$', fontsize=20)
    ax.set_xlabel(r'$\frac{\alpha_{inner}}{\alpha_{outer}}$', fontsize=20)
    if over_1:
        ax.scatter([2],[3],c='r')
def ratios_analysis(df,plot=True):
    a = np.arange(1,10,0.1)
    b = np.arange(1,10,0.1)
    A, B = np.meshgrid(a, b, indexing='ij')
    df_tsp = pd.read_csv(name)
    shape = np.shape(A)
    Z = np.zeros_like(A)
    a_out = 1
    b_out = 1
    for i in range(shape[0]):
        for j in range(shape[1]):
            a_in = A[i,j]/a_out
            b_in = B[i,j]/b_out
            df_new = re_calculate_complexity(df,a_in,a_out,b_in,b_out)
            m = analysis(df_new)
            Z[i,j] = m['coef']['b']
    if plot:
        # plt.rcParams['text.usetex'] = True
        plot_contour_ratio(A,B,Z)
    return (A,B,Z)
def ratios_analysis_under_1(df,plot=True):
    a = np.arange(0.01,1,0.01)
    b = np.arange(0.01,1,0.01)
    A, B = np.meshgrid(a, b, indexing='ij')
    df_tsp = pd.read_csv(name)
    shape = np.shape(A)
    Z = np.zeros_like(A)
    a_out = 1
    b_out = 1
    for i in range(shape[0]):
        for j in range(shape[1]):
            a_in = A[i,j]/a_out
            b_in = B[i,j]/b_out
            df_new = re_calculate_complexity(df,a_in,a_out,b_in,b_out)
            m = analysis(df_new)
            Z[i,j] = m['coef']['b']
    if plot:
        # plt.rcParams['text.usetex'] = True
        plot_contour_ratio(A,B,Z)
    return (A,B,Z)            

df = load_results()
# Alpha,Beta,Z = ratios_analysis(df)
# Alpha_u,Beta_u,Z_u = ratios_analysis_under_1(df)
# plot_contour_ratio(Alpha,Beta,Z)
# results = analysis_with_plot(df)


# df_new = re_calculate_complexity(df,.2,.1,.3,.1)
# results = analysis_with_plot(df_new)
# df_tsp = pd.read_csv(name)
# out = []
# for i in range(len(df_tsp)):
#     inner = df_tsp['interior points'].iloc[i]
#     outer = df_tsp['outside points'].iloc[i]
#     out.append(complexity_tsp_given_points(inner,outer,2,1,2,1))
# df_tsp['New_C'] = out
# df['C']= df['id'].map(df_tsp[['id','New_C']].set_index('id').squeeze())
# m = analysis(df)

class TSP_Analysis(Monte_Carlo):
    def __init__(self):
        self.x_label = 'exponent coeficient'
        self.df = load_results()
        self.df_tsp = pd.read_csv(name)
        self.alpha_inner = variable(triangular,0,3,2) 
        self.alpha_outer = variable(triangular,0,2,1) 
        self.beta_inner = variable(triangular,0,3,2)
        self.beta_outer = variable(triangular,0,2,1)
        self.deterministic = analysis(self.df)['coef']['b']
        self.history = []
    def box_plot_Complexity(self):
        fig = plt.figure()
        hist = self.history
        C = []
        for i in hist:
            C.append(i['Problem Complexity'])
        cols = get_ids()
        cols[30] = 'id 30'
        C = pd.DataFrame(C)
        C =C.drop(columns=[29])
        C = C.rename(columns=cols)
        bp = C.boxplot()
        plt.xlabel('Problem ID')
        plt.ylabel('Complexity')
        return C
    def plot_sensitivity(self):
        df = self.sensitivity_df
        idx = df.index.to_list()
        n = len(idx)
        pos = np.arange(n) + .5 # bars centered on the y axis
        fig, (ax_left, ax_right) = plt.subplots(ncols=2)
        ax_left.barh(pos, df[1]*100, align='center', facecolor='cornflowerblue')
        ax_right.set_yticks([])
        ax_left.spines['top'].set_visible(False)
        ax_left.spines['right'].set_visible(False)
        ax_left.spines['bottom'].set_visible(False)
        ax_left.spines['left'].set_visible(False)
        ax_left.xaxis.set_major_formatter(mtick.PercentFormatter())
        # ax_left.set_xlabel('Complexity change')
        # ax_left.invert_xaxis()
        ax_right.barh(pos, df[0]*100, align='center', facecolor='limegreen')
        ax_left.set_yticks(pos)
        idx_map = {'alpha_in':r'$\alpha_{inner}$','alpha_out':r'$\alpha_{outer}$',
                   'beta_in':r'$\beta_{inner}$','beta_out':r'$\beta_{outer}$'}
        idx_label = map(lambda x: idx_map[x], idx)
        ax_left.set_yticklabels(idx_label, ha='center', x=-0.05,fontsize=15)
        
        ax_right.spines['top'].set_visible(False)
        ax_right.spines['right'].set_visible(False)
        ax_right.spines['bottom'].set_visible(False)
        ax_right.spines['left'].set_visible(False)
        ax_right.xaxis.set_major_formatter(mtick.PercentFormatter())
        fig.subplots_adjust(wspace=0.0)
        fig.suptitle('Complexity Sensitivity')
        fig.show()
    def sensitivity(self):
        det_alpha_inner = 2
        det_alpha_outer = 1 
        det_beta_inner = 3
        det_beta_outer = 1
        delta = 0.10
        alpha_in_plus = []
        alpha_in_minus =[]
        out_normal = []
        alpha_out_plus = []
        alpha_out_minus =[]
        beta_in_plus = []
        beta_in_minus =[]
        beta_out_plus = []
        beta_out_minus =[]
        for i in range(len(self.df_tsp)):
            inner = self.df_tsp['interior points'].iloc[i]
            outer = self.df_tsp['outside points'].iloc[i]
            alpha_in_plus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner*(1+delta),det_alpha_outer,
                                                   det_beta_inner,det_beta_outer))
            alpha_in_minus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner*(1-delta),det_alpha_outer,
                                                   det_beta_inner,det_beta_outer))
            out_normal.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner,det_alpha_outer,
                                                   det_beta_inner,det_beta_outer))
            
            alpha_out_plus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner,det_alpha_outer*(1+delta),
                                                   det_beta_inner,det_beta_outer))
            alpha_out_minus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner,det_alpha_outer*(1-delta),
                                                   det_beta_inner,det_beta_outer))
            
            beta_in_plus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner,det_alpha_outer,
                                                   det_beta_inner*(1+delta),det_beta_outer))
            beta_in_minus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner,det_alpha_outer,
                                                   det_beta_inner*(1-delta),det_beta_outer))
            
            beta_out_plus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner,det_alpha_outer,
                                                   det_beta_inner,det_beta_outer*(1+delta)))
            beta_out_minus.append(complexity_tsp_given_points_full(inner,outer,
                                                   det_alpha_inner,det_alpha_outer,
                                                   det_beta_inner,det_beta_outer*(1-delta)))
        normal = pd.DataFrame(out_normal)
        alpha_in_plus = pd.DataFrame(alpha_in_plus)
        alpha_in_minus = pd.DataFrame(alpha_in_minus)
        alpha_out_plus = pd.DataFrame(alpha_out_plus)
        alpha_out_minus = pd.DataFrame(alpha_out_minus)
        
        beta_in_plus = pd.DataFrame(beta_in_plus)
        beta_in_minus = pd.DataFrame(beta_in_minus)
        beta_out_plus = pd.DataFrame(beta_out_plus)
        beta_out_minus = pd.DataFrame(beta_out_minus)
        
        alpha_in_plus = (alpha_in_plus-normal)/normal
        alpha_in_minus = (alpha_in_minus-normal)/normal
        alpha_out_plus = (alpha_out_plus-normal)/normal
        alpha_out_minus = (alpha_out_minus-normal)/normal
        
        beta_in_plus = (beta_in_plus-normal)/normal
        beta_in_minus = (beta_in_minus-normal)/normal
        beta_out_plus = (beta_out_plus-normal)/normal
        beta_out_minus = (beta_out_minus-normal)/normal
        
        out = pd.DataFrame({'beta_in':[beta_in_plus.mean()['C'],beta_in_minus.mean()['C']],
               'beta_out':[beta_out_plus.mean()['C'],beta_out_minus.mean()['C']],
               'alpha_in':[alpha_in_plus.mean()['C'],alpha_in_minus.mean()['C']],
               'alpha_out':[alpha_out_plus.mean()['C'],alpha_out_minus.mean()['C']]}).T
        out = out.sort_values(by=0,ascending = True)
        self.sensitivity_df = out
        self.plot_sensitivity()
        return {'normal':normal,'alpha_in_plus':alpha_in_plus,'alpha_in_minus':alpha_in_minus,'alpha_out_plus':alpha_out_plus,'alpha_out_minus':alpha_out_minus,
                'beta_in_plus':beta_in_plus,'beta_in_minus':beta_in_minus,'beta_out_plus':beta_out_plus,'beta_out_minus':beta_out_minus,'out':out}
    
        
        
    def evaluate(self):
        out = []
        a_in = self.alpha_inner()
        a_out = self.alpha_outer()
        b_in = self.beta_inner()
        b_out = self.beta_outer()
        for i in range(len(self.df_tsp)):
            inner = self.df_tsp['interior points'].iloc[i]
            outer = self.df_tsp['outside points'].iloc[i]
            out.append(complexity_tsp_given_points(inner,outer,
                                                   a_in,a_out,
                                                   b_in,b_out))
        self.df_tsp['New_C'] = out
        self.df['C']= self.df['id'].map(self.df_tsp[['id','New_C']].set_index('id').squeeze())
        
        m = analysis(self.df)
        m['Complexity Coef']={}
        m['Complexity Coef']['alpha_inner']= a_in
        m['Complexity Coef']['alpha_outer']= a_out
        m['Complexity Coef']['beta_inner']= b_in
        m['Complexity Coef']['beta_outer']= b_out
        m['Problem Complexity'] = out
        self.history.append(m)
        return m['coef']['b']

# tsp = TSP_Analysis()
# out = tsp.sensitivity()
# df_p = pd.DataFrame(out['plus'])
# df_m = pd.DataFrame(out['minus'])
# df_n = pd.DataFrame(out['normal'])
# results = tsp.simulation(10000)
# tsp.plot_cdf()
# hist = tsp.history
# bp = tsp.box_plot_Complexity()
# df = load_results()
# df['tn']=df["time"] * df["normalized_cost"]

# tf = np.invert(np.isnan(df['time']))
# df2 = df[tf].groupby('C').mean()

# df2 = df[['C','id']][tf].groupby('id').mean()
# df2['count'] = df[['C','id']][tf].groupby('id').count()
# df2['time'] = df[['id','time']][tf].groupby('id').mean()
# df2['timeSTD'] = df[['id','time']][tf].groupby('id').std()
# df2['cost'] = df[['id','cost']][tf].groupby('id').mean()
# df2['costSTD'] = df[['id','cost']][tf].groupby('id').std()
# df2['normalized_cost'] = df[['id','normalized_cost']][tf].groupby('id').mean()
# df2['normalized_costSTD'] = df[['id','normalized_cost']][tf].groupby('id').std()
# df2['tn'] = df[['id','tn']][tf].groupby('id').mean()
# df2['tnSTD'] = df[['id','tn']][tf].groupby('id').std()
# df3=df2.copy()
# for i in range(len(df2)):
#     cur = df2.iloc[i]
#     count = cur['count']
#     for j in range(int(count)-1):
#         df3= df3.append(cur) 
# x = df2.index
# y = df2['tn']
# metric = Custom_Fit(x, y, fit_func_2)