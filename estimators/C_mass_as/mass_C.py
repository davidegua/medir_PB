import numpy as np
import ctypes
from ctypes import *
from numpy.ctypeslib import ndpointer

"define a pointer for 1D arrays"
_doublep  = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')
_longp    = ndpointer(ctypes.c_long, flags='C_CONTIGUOUS')
"define a pointer for 1D arrays INT "
_intp  = ndpointer(ctypes.c_int, flags='C_CONTIGUOUS')
"define a pointer for 2D arrays"
_doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')
"function to convert 2D array into a proper format for C"
def c_2d_inp(x):
    return (x.__array_interface__['data'][0]
            + np.arange(x.shape[0]) * x.strides[0]).astype(np.uintp)

path_lib = './estimators/C_mass_as/'
dftlib  = CDLL(path_lib + 'mass_as.so')

dftlib.mass_ass_pq5s.restype  = None
dftlib.mass_ass_pq5s.argtypes = [_doublep, _doublep, _doublep,_doublep,_doublep, c_double, c_double, c_int,  c_int]

def mass_ass_pq5s(xyzw,L2,L1,Ngrid):

    delta = np.zeros(Ngrid**3,dtype=float)

    dftlib.mass_ass_pq5s(delta, np.array(xyzw[:,0]), np.array(xyzw[:,1]), np.array(xyzw[:,2]), np.array(xyzw[:,3]),
                          L2,L1,Ngrid,int(xyzw.shape[0]))

    return delta