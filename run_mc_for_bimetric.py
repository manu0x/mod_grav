#!/usr/bin/env python
# coding: utf-8




import ctypes
import sys
import os 
import numpy as np
from ctypes import POINTER
import matplotlib.pyplot as plt
from mc_for_bimetric import log_likelihood,log_post
import emcee



cwd =  os.getcwd()
lib_h = ctypes.CDLL(cwd + "/friedmann_Cs_lib.so") 




def load_data(fname="data.txt"):
    f = np.genfromtxt(fname,delimiter="\t")
    z = np.squeeze(f[:,0])
    fs8_data = np.squeeze(f[:,1])
    fs8_sig2 = np.square(np.squeeze(f[:,2]))
    om_fid = np.squeeze(f[:,3])
    
    c = np.diag(fs8_sig2)
    c[5,5] = 6.4000*0.001
    c[4,4] = 3.969*0.001
    c[3,3] = 5.184*0.001
    c[4,5] = 2.570*0.001
    c[3,5] = 0.000
    c[5,4] = 2.570*0.001
    c[3,4] = 2.540*0.001
    c[5,3] = 0.000*0.001
    c[4,3] = 2.540*0.001
    
    c_i = np.linalg.inv(c)
    rat_den = np.sqrt(om_fid*np.power(1.0+z,3.0)+(1.0-om_fid))
    mz = -1.0*z
    sa = mz.argsort()
    sa_inv = mz.argsort().argsort()
    
    cn = np.array([1,1,2,1,2,1,1,1,1,1,1,1,1,1,1,3,3,3])
    
    return z,fs8_data,c_i,rat_den,sa,sa_inv,cn
    
    



z,f,ci,r,s,si,cn = load_data()





l = log_post([0.4,0.0,0.79],z,f,r,ci,cn,s,si)




nwalkers = 100
ndim = 3
ini_guess =  np.array([0.3,0.0,0.79])+1e-4 * np.random.randn(nwalkers, ndim)


filename = "run.h5"
backend = emcee.backends.HDFBackend(filename)
backend.reset(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_post, args=(z,f,r,ci,cn,s,si))



sampler.run_mcmc(ini_guess, 100, progress=True);




samples =sampler.get_chain()





