import numpy as np
import pandas as pd
from estimators.fftw import pydft
from estimators.C_mass_as.mass_C import mass_ass_pq5s
import estimators.pyliansfun.readgadget as readgadget
#from astropy.io import ascii
import random

"""
save measurements
"""
def save_measu(out_pk,out_bk,Ngrid,Lbox,sdetails,sim_idx,snap,dkfac,Nlacing,measure_bk,fraction):
    "save measurements"
    """POWER SPECTRUM"""
    guadatapk = np.array([out_pk['kef'],out_pk['p0k'],out_pk['p0k_msn'],out_pk['kcb']]).T

    np.savetxt('./output/guaspec_pkout_%s_%dNgrid_%dLbox_%drel_%dsnap_%ddkf_%dNlac_%dfractra.txt'
               % (sdetails, Ngrid, Lbox, sim_idx, snap, int(dkfac * 100), Nlacing, int(fraction * 100)),
               guadatapk,fmt='%.5e')


    """BISPECTRUM"""
    if measure_bk==True:
        guadatabk = np.array([out_bk['xbk'], out_bk['b123'], out_bk['b123_msn'],
                                  out_bk['tr_co_eff'][:, 0], out_bk['tr_co_eff'][:, 1], out_bk['tr_co_eff'][:, 2]]).T

        np.savetxt('./output/guaspec_bkout_%s_%dNgrid_%dLbox_%drel_%dsnap_%ddkf_%dNlac_%dfractra.txt'
                   % (sdetails, Ngrid, Lbox, sim_idx, snap, int(dkfac * 100), Nlacing, int(fraction * 100)),
                      guadatabk, fmt='%.5e')


"""
indexing of the box for the Fourier transform
"""
def box_ordering(Ngrid):
    dimtot = Ngrid ** 3
    lspace = np.linspace(0,dimtot-1,dimtot,dtype=int)

    iar =  lspace // Ngrid**2
    jar = (lspace %  Ngrid**2) // Ngrid
    kar = (lspace %  Ngrid**2) %  Ngrid

    Iar = np.copy(iar)
    Jar = np.copy(jar)
    Kar = np.copy(kar)

    Iar[iar>(Ngrid/2)] = - (Ngrid - iar[iar>(Ngrid/2)])
    Jar[jar>(Ngrid/2)] = - (Ngrid - jar[jar>(Ngrid/2)])
    Kar[kar>(Ngrid/2)] = - (Ngrid - kar[kar>(Ngrid/2)])

    Lar   = Iar**2 + Jar**2 + Kar**2
    Lared = Lar[Lar<=(Ngrid**2/4.)]
    lared = lspace[Lar<=(Ngrid**2/4.)]


    orderedL = np.argsort(Lared)
    lLar = np.vstack((lared[orderedL],Lared[orderedL])).T

    del iar,jar,kar,Iar,Jar,Kar,Lar,Lared,lared,orderedL
    return lLar


"""""""""""""""""""""""""""""""""""""""""""""""
START FROM THE PARTICLES
"""""""""""""""""""""""""""""""""""""""""""""""
"""
many realisation load data function for periodic boxes
"""
def load_data_and_cstep_manyrel( Ngrid, deltak_fac,snap,path,Nlacing,sdetails,fraction):

    name_file = path + "snap_00%d"%snap
    xyzw, Lbox = Qui_gadg(name_file,False,fraction)

    print("Lbox %d Mpc/h"%Lbox)
    kf = 2. * np.pi / Lbox
    Lmin = 0.
    Lmax = Lbox

    deltak = kf * deltak_fac
    kmaxpk = kf * (Ngrid / 2. - deltak_fac)
    kmaxbk = kf * (Ngrid / 3. - deltak_fac)
    kmaxtk = kf * (Ngrid / 6. - deltak_fac)  # not four but half of the kmaxbk since for the SN we need bk for k=2*kmaxtk

    "arbitrary to reduce memory load"
    kmaxpk = 0.12
    kmaxbk = 0.12
    kmaxtk = 0.06

    Ndata = np.float64(xyzw.shape[0])

    ngperiodic = Ndata/Lbox**3.
    print("average number density of galaxies periodic box %.5f"%ngperiodic)

    print ("shape input data ", xyzw.shape)

    kx = np.zeros(Ngrid)
    ir = np.arange(Ngrid)

    ig = np.arange(Ngrid**3/2 + Ngrid**2)

    i1 = ir < (Ngrid/2 + 1)
    i2 = ir > Ngrid/2
    kx[i1] = ir[i1] * kf
    kx[i2] = - (Ngrid - ir[i2]) * kf

    for i in range(1,Nlacing+1):
        print("interlacing step %d"%i)

        L2 = Lmax - (Lmax - Lmin)/Ngrid*1.*1./Nlacing*1.*(i-1)
        L1 = Lmin - (Lmax - Lmin)/Ngrid*1.*1./Nlacing*1.*(i-1)

        deltax = mass_ass_pq5s(xyzw,L2,L1,Ngrid)

        deltax /= Ndata
        deltax -=  1. / np.float64(Ngrid) ** 3

        print("before performing fft")
        "perform FFT"
        if i == 1:
            dkre, dkim  = pydft.dx_to_dk(deltax, Ngrid, L2 - L1, mode_mass_ass=6)

            ii = np.array(ig / (Ngrid ** 2. / 2. + Ngrid), dtype=np.int32)
            jj = np.array((ig - ii * (Ngrid ** 2 / 2. + Ngrid)) / (Ngrid / 2. + 1.), dtype=np.int32)
            kk = np.array(ig - ii * (Ngrid ** 2 / 2. + Ngrid) - jj * (Ngrid / 2. + 1.), dtype=np.int32)

        else:
            dkreb,dkimb = pydft.dx_to_dk(deltax, Ngrid, L2 - L1, mode_mass_ass=6)

            phase_cos=np.cos(((Lmax-Lmin)/np.float64(Ngrid)) * (i * 1.-1.)/np.float64(Nlacing) * (kx[ii] + kx[jj] + kx[kk]))
            phase_sin=np.sin(((Lmax-Lmin)/np.float64(Ngrid)) * (i * 1.-1.)/np.float64(Nlacing) * (kx[ii] + kx[jj] + kx[kk]))

            new_deltak_re_b = dkreb * phase_cos - dkimb * phase_sin
            new_deltak_im_b = dkreb * phase_sin + dkimb * phase_cos

            dkre += new_deltak_re_b
            dkim += new_deltak_im_b

    dkre /= np.float64(Nlacing)
    dkim /= np.float64(Nlacing)

    lLar = box_ordering(Ngrid)

    cstep_d, Ik_x, Jk_x = common_step(dkre,dkim,Lbox,Ngrid,deltak,kmaxpk,lLar,ngperiodic)

    cstep_d['kmaxbk']  = kmaxbk
    cstep_d['kmaxtk']  = kmaxtk
    cstep_d['s_deta']  = sdetails

    print("real space")
    return cstep_d, Ik_x, Jk_x


"""
prepare a common box with the k-indexes binned in the same way
start from the catalogue, read particles and compute the delta_x,
delta_k, and finally Ikx, Jkx fields
"""
def common_step(dkre,dkim, Lbox, Ngrid,deltak,kmaxpk,lLar,ng):

    I22 = Lbox ** (-3.)
    I33 = Lbox ** (-6.)
    I44 = Lbox ** (-9.)

    kf = 2*np.pi / np.float64(Lbox)

    kny_pkbktk = kf * Ngrid * np.array([1./2.,1./3.,1./4.])
    print("Ny. freq pk %.3f bk %.3f tk %.3f"%(kny_pkbktk[0],kny_pkbktk[1],kny_pkbktk[2]))

    Nbins = np.floor(kny_pkbktk[0]/deltak).astype(int)
    Nmax  = int(np.floor(kmaxpk/deltak) + 1)

    if (Nmax>Nbins): print("Nbins %d smaller than Nmax %d !!!!"%(Nbins,Nmax))
    print("given the chosen kmax for pk %.5f and bin size %.5f , there are %d bins"%(kmaxpk,deltak,Nmax))

    dim_r2c = np.int64(Ngrid**3/2 + Ngrid**2)

    print("--- calculating I_k(x) shells ---")
    Ik_x, Ikeff, Ineff = pydft.c2r_dk_to_Ikx_v2(dkre, dkim, lLar, Ngrid, Nmax, deltak, kf)
    print("--- calculating J_k(x) shells ---")
    Jk_x, Jkeff, Jneff = pydft.c2r_dk_to_Ikx_v2(np.ones(dim_r2c, dtype=np.float64),
                                                np.zeros(dim_r2c,dtype=np.float64),
                                                lLar, Ngrid, Nmax, deltak, kf)

    print("\n Ikeff - Keff \n")
    print(Ikeff)

    print("creating dictionary")
    "create_dictionary"
    cstep_d = {}
    cstep_d['deltak'] = deltak
    cstep_d['Ngrid']  = Ngrid
    cstep_d['I22'] = I22
    cstep_d['I33'] = I33
    cstep_d['I44'] = I44
    cstep_d['k_f'] = kf
    cstep_d['kny_pkbktk'] = kny_pkbktk
    cstep_d['Nmax'] = Nmax
    cstep_d['keff'] = Ikeff
    cstep_d['ng']   = ng
    cstep_d['Lbox'] = Lbox

    return cstep_d, Ik_x, Jk_x


"function to read gadget files from Quijote simulations"
#adapted from Quijote's documentation, francisco villaescusa's github
def Qui_gadg(filepath,rsd,fraction):
    # input files
    ptype    = [1] #[1](CDM), [2](neutrinos) or [1,2](CDM+neutrinos)

    print("filepath: ",filepath)
    # read header
    header   = readgadget.header(filepath)
    BoxSize  = header.boxsize/1e3  #Mpc/h
    Nall     = header.nall         #Total number of particles
    Masses   = header.massarr*1e10 #Masses of the particles in Msun/h
    Omega_m  = header.omega_m      #value of Omega_m
    Omega_l  = header.omega_l      #value of Omega_l
    h        = header.hubble       #value of h
    redshift = header.redshift     #redshift of the snapshot
    Hubble   = 100.0*np.sqrt(Omega_m*(1.0+redshift)**3+Omega_l)#Value of H(z) in km/s/(Mpc/h)

    print("Lbox = %.0f Mpc/h"%BoxSize)

    # read positions, velocities and IDs of the particles
    pos = readgadget.read_block(filepath, "POS ", ptype)/1e3  #positions in Mpc/h
    vel = readgadget.read_block(filepath, "VEL ", ptype)     #peculiar velocities in km/s
    ids = readgadget.read_block(filepath, "ID  ", ptype)-1   #IDs starting from 0

    print("object read : %d"%pos.shape[0])

    if rsd == True:
        scafa = 1./(1.+ redshift)
        print("converting to redshift space ...")

        pos[:,2] += vel[:,2]/(Hubble * scafa)

        pos[pos[:,2]>BoxSize,2] = (pos[pos[:,2]>BoxSize,2] + 1.*BoxSize)%(1.*BoxSize)
        pos[pos[:,2]<0.0,2] = (pos[pos[:,2]<0.0,2] + 1.*BoxSize)%(1.*BoxSize)

    "shot noise test sub-select set of particles"
    Nsub = np.int32(fraction * pos.shape[0])
    idx  = np.arange(pos.shape[0])
    random.shuffle(idx)
    pos  = pos[idx[:Nsub]]

    weigh = np.ones([pos.shape[0],1])
    xyzw  = np.hstack((pos,weigh))

    print(r"Number of objects:  %d, %d per cent  of the total"%(xyzw.shape[0],fraction*100))

    return xyzw, BoxSize














