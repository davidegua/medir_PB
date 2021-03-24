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

path_lib = './estimators/fftw/'
dftlib  = CDLL(path_lib + 'dft.so')

dftlib.r2c_dx_to_dk.restype  = None
dftlib.r2c_dx_to_dk.argtypes = [_doublep, _doublep, _doublep, c_int, c_double, c_double, c_int]

def dx_to_dk(deltax,ngrid,L,mode_mass_ass):

    dx1d = deltax.ravel()

    compdim = int(0.5 * ngrid**3 + ngrid**2)
    dkre = np.zeros(compdim,dtype=np.float64)
    dkim = np.zeros(compdim,dtype=np.float64)

    dftlib.r2c_dx_to_dk(dx1d,dkre,dkim,ngrid,0.,L,mode_mass_ass)

    return dkre,dkim


# "quadrupole and hexadecapole delta k field using z-axis as LOS"
# dftlib.r2c_dx_to_dk_multip.restype  = None
# dftlib.r2c_dx_to_dk_multip.argtypes = [_doublep, _doublep, _doublep, _doublep, _doublep, c_int, c_double, c_double, c_int]
#
# def dx_to_dk_multip(deltax,ngrid,L,mode_mass_ass):
#
#     dx1d = deltax.ravel()
#
#     compdim = int(0.5 * ngrid**3 + ngrid**2)
#     dkreq = np.zeros(compdim,dtype=np.float64)
#     dkimq = np.zeros(compdim,dtype=np.float64)
#     dkreh = np.zeros(compdim, dtype=np.float64)
#     dkimh = np.zeros(compdim, dtype=np.float64)
#
#     dftlib.r2c_dx_to_dk(dx1d,dkreq,dkimq,dkreh,dkimh,ngrid,0.,L,mode_mass_ass)
#
#     return dkreq,dkimq,dkreh,dkimh
#

"version same as Rustico algorithm"
dftlib.c2r_dk_to_Ikx_v2.restype  = None
dftlib.c2r_dk_to_Ikx_v2.argtypes = [_doublep, _doublep, _doublepp, _longp, _longp, #_longp,
                                    c_int, c_int,c_double, c_double, _doublep, _intp]

# def c2r_dk_to_Ikx_v2(deltak_re,deltak_im,lLar,idL,ngrid,nbmax,Deltak,kf):
def c2r_dk_to_Ikx_v2(deltak_re,deltak_im,lLar,ngrid,nbmax,Deltak,kf):
    realdim = int(ngrid**3)
    Neffar  = np.zeros(nbmax,dtype=np.int32)
    Keffar  = np.zeros(nbmax,dtype=np.float64)
    Ik_x    = np.zeros([nbmax,realdim],dtype=np.float64)

    print("before c function")
    dftlib.c2r_dk_to_Ikx_v2(deltak_re,deltak_im,c_2d_inp(Ik_x),
                            np.array(lLar[:,0],dtype=np.int64),np.array(lLar[:,1],dtype=np.int64),#np.array(idL,dtype=long),
                            ngrid,nbmax,Deltak,kf,Keffar,Neffar)
    print("outside c function")
    return Ik_x, Keffar, Neffar



