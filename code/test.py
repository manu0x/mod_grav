

import os
import numpy as np

code_dir  = (os.path.dirname(os.path.realpath(__file__)))

#cwd =  os.getcwd()
#lib_h = ctypes.CDLL(cwd + "/friedmann_Cs_lib.so") 



print(data_dir,code_dir,fname) 


f = np.genfromtxt(fname,delimiter="\t")
print(f.shape)
