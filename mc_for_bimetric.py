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
[ctypes.c_double,ctypes.c_double,ctypes.c_double,\
POINTER(ctypes.c_double),POINTER(ctypes.c_double),POINTER(ctypes.c_double),POINTER(ctypes.c_long),ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_double]
lib_h.restype = ctypes.c_int





def log_likelihood(theta,z,data,rat_den,C_I,cn,s,si):
   
  
    omega_dm_0 = theta[0]
    B1 = theta[1]
    s8 = theta[2]
    n = z.shape[0]
    zz = ctypes.c_double*(n)
    ratio_den = ctypes.c_double*(n)
    cnw = ctypes.c_int*(n)
    fs8 = ctypes.c_double*n

    zw = z[s]
    rat_denw = rat_den[s]
    cnww = cn[s]
    
    cnw = np.ctypeslib.as_ctypes(cnww) 
    
    zz = np.ctypeslib.as_ctypes(zw) 
    ratio_den = np.ctypeslib.as_ctypes(rat_denw) 
    fs8 = np.ctypeslib.as_ctypes(np.zeros((n,)))
    #print(cnw[1:10],cn)

    l = lib_h.bimetric_fs8_log_like(omega_dm_0,B1,s8,zz,fs8,ratio_den,cnw,n,0,1,1.0)  
    fs8 = np.ctypeslib.as_array(fs8)
    data = np.ctypeslib.as_array(data)
    fs8 = np.reshape(fs8,(fs8.shape[0],1))
    data = np.reshape(data,(data.shape[0],1))
    fs8 = fs8[si]
    diff = data-fs8
    
    lklhd = -np.matmul(np.matmul(np.transpose(diff),C_I),diff)
    
    return lklhd[0,0]





def log_prior(theta):
    omega_dm_0 = theta[0]
    B1 = theta[1]
    s8 = theta[2]

    if 0.01<omega_dm_0<0.9 and 0.0<=B1<=6.0 and 0.01<s8<1.0:
        return 0.0
    else:
        return -np.inf


def log_post(theta,z,data,rat_den,C_I,cn,s,si):
    log_pr = log_prior(theta)
    if not(np.isfinite(log_pr)):
        return -np.inf
    else:
        return log_pr + log_likelihood(theta,z,data,rat_den,C_I,cn,s,si)

