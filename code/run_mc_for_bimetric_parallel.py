#!/usr/bin/env python
# coding: utf-8




import ctypes
import sys
import os 
import time
import numpy as np
from ctypes import POINTER
import matplotlib.pyplot as plt
from mc_for_bimetric import log_likelihood,log_post
import emcee
from schwimmbad import MPIPool
from mpi4py import MPI



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





#l = log_post([0.4,0.0,0.79],z,f,r,ci,cn,s,si)






ndim = 3
n_steps = int(sys.argv[1])
nwalkers = int(sys.argv[2])
print(n_steps,nwalkers)






with MPIPool() as pool:
	if not pool.is_master():
		pool.wait()
		sys.exit(0)

	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()

	print(rank)
	
	np.random.seed(rank)

	filename = "run.h5"
	backend = emcee.backends.HDFBackend(filename)
	backend.reset(nwalkers, ndim)
	


	
	ini_guess =  np.array([0.3,0.0,0.79])+1e-4 * np.random.randn(nwalkers, ndim)


	

	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_post, pool=pool,backend = backend, args=(z,f,r,ci,cn,s,si))

	start = time.time()

	sampler.run_mcmc(ini_guess, n_steps, progress=True);

	end = time.time()
	
	print(end - start)










