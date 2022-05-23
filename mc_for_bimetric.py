#!/usr/bin/env python
# coding: utf-8




import ctypes
import sys
import os 
import numpy as np
from ctypes import POINTER
import matplotlib.pyplot as plt


# In[2]:


cwd =  os.getcwd()





lib_h = ctypes.CDLL(cwd + "/friedmann_Cs_lib.so") 





lib_h.bimetric_fs8_log_like.argtypes=\
[ctypes.c_double,ctypes.c_double,ctypes.c_double,POINTER(ctypes.c_double),POINTER(ctypes.c_double),ctypes.c_int,ctypes.c_int,ctypes.c_int]
lib_h.restype = ctypes.c_int





def calc_likelihood(z,data,C_I,omega_dm_0=0.3,B1=0.0,s8=0.79):
    n = z.shape[0]
    zz = ctypes.c_double*(n)
    zz = np.ctypeslib.as_ctypes(z) 
    fs8 = ctypes.c_double*n
    fs8 = np.ctypeslib.as_ctypes(np.zeros((n,)))
    l = lib_h.bimetric_fs8_log_like(omega_dm_0,B1,s8,zz,fs8,n,1,1)
    fs8 = np.reshape(fs8,(fs8.shape[0],1))
    data = np.reshape(data,(data.shape[0],1))
    diff = data-fs8
    print(fs8.shape,diff.shape)
    lklhd = -np.matmul(np.matmul(np.transpose(diff),C_I),diff)
    
    return lklhd





