# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:04:37 2021

@author: Ricardo Hopker
"""

import numpy as np
import xlwings as xw
from complexFuncs import complexity, DSM_rearrange,DSM_idx_arrange,DSM_colapse_simple,DSM_colapse_complex,DSM_colapse,DSM_colapse2

import pandas as pd

@xw.func
@xw.arg('DSM',np.array)
@xw.ret(expand='table')
def complexityXL(DSM):
    # print(DSM)
    out = complexity(DSM)
    # print(out)
    return out['C']
@xw.func
@xw.arg('DSM',np.array)
@xw.ret(expand='table')
def complexityXL_Full(DSM):
    # print(DSM)
    out = complexity(DSM)
    # print(out)
    outArr = [
        ['C',out['C']],
             ['C1',out['C1']],
             ['C2',out['C2']],
             ['C3',out['C3']]]
    return outArr

@xw.func
@xw.arg('DSM',np.array)
@xw.arg('cols',np.array,dtype='int')
@xw.ret(expand='table')
def DSM_rearrangeXL(DSM,cols):
    # print(DSM)
    out = DSM_rearrange(DSM,cols)
    # print(out)
    return out

@xw.func
@xw.arg('DSM',np.array)
@xw.arg('idx',np.array,dtype='int')
@xw.ret(expand='table')
def DSM_idx_arrange_XL(DSM,idx):
    out = DSM_idx_arrange(DSM,idx)
    # print(out)
    return out

@xw.func
@xw.arg('DSM',np.array)
@xw.arg('idx',np.array,dtype='int')
@xw.ret(expand='table')
def DSM_colapse_complex_XL(DSM,idx):
    out = DSM_colapse_complex(DSM,idx)
    # print(out)
    return out

@xw.func
@xw.arg('DSM',np.array)
@xw.arg('args',np.array,dtype='int')
@xw.ret(expand='table')
def DSM_multiple_colapse_complex_XL(DSM,*args):
    idx = []
    for x in args:
        idx.append(x)
    # print(idx)
    out = DSM_colapse2(DSM,idx)
    # print(out)
    return out

@xw.func
@xw.arg('DSM',np.array)
@xw.arg('idx',np.array,dtype='int')
@xw.ret(expand='table')
def DSM_colapse_simple_XL(DSM,idx):
    out = DSM_colapse_simple(DSM,idx)
    # print(out)
    return out
@xw.func
@xw.arg('DSM',np.array)
@xw.arg('args',np.array,dtype='int')
@xw.ret(expand='table')
def DSM_multiple_colapse_simple_XL(DSM,*args):
    idx = []
    for x in args:
        idx.append(x)
    # print(idx)
    out = DSM_colapse(DSM,idx)
    # print(out)
    return out
if __name__ == '__main__':
    xw.serve()