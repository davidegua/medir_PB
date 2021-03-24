import numpy as np

"""
 measure the power spectrum
"""
def Pk_periodic(info_d,Ikx,Jkx):

    kcb     = np.zeros(info_d['Nmax'])
    p0k     = np.zeros(info_d['Nmax'])
    p0k_msn = np.zeros(info_d['Nmax'])
    Nk      = np.zeros(info_d['Nmax'])
    ng      = info_d['ng']

    for j in range(info_d['Nmax']):
        if j%5 ==0: print(j)

        kcb[j] = (j+1) * info_d['deltak']
        "normalisation"
        Nk[j] = np.einsum('i,i',Jkx[j],Jkx[j]) / info_d['Ngrid'] ** 3.

        "Ikx product"
        p0k[j] = np.einsum('i,i',Ikx[j],Ikx[j]) / info_d['Ngrid'] ** 3.
        p0k[j] = p0k[j] /float(Nk[j]) /info_d['I22']

        p0k_msn[j] = p0k[j] - 1./ng

    out_pk = {}
    out_pk['kcb']     = kcb
    out_pk['kef']     = info_d['keff']
    out_pk['p0k']     = p0k
    out_pk['p0k_msn'] = p0k_msn
    out_pk['Nk']      = Nk

    return out_pk


"""
measure the bispectrum
"""
def Bk_periodic(info_d,kmin,kmax,p0k_msn,Ikx,Jkx):
    ng     = info_d['ng']
    Nmaxbk = int((kmax-kmin)/info_d['deltak']) - 1

    print("Nmaxbk %d"%Nmaxbk, "kmin %.3f"%kmin)

    counts = np.zeros((Nmaxbk, Nmaxbk, Nmaxbk), dtype=float) #default double prec
    bkijl  = np.zeros((Nmaxbk, Nmaxbk, Nmaxbk), dtype=float)

    triangles, indexes_tr, triangles_eff = [], [],[]
    p0k_i, p0k_j, p0k_l   = [], [], []
    b123_arr, b123_msn_arr, q123_arr, bksn_arr = [], [], [], []


    tr_config_file_name = "./configurations/bk_%dNgrid_%ddkfac_%dkmax_%dLbox.npz"\
                          %(info_d['Ngrid'],int(info_d['deltak']/info_d['k_f']),int(kmax*100),int(info_d['Lbox']))
    try :
        indexes_tr = np.load(tr_config_file_name)['indexes_tr']
        print("found wisdom file with tr coords for %d triangles!"%indexes_tr.shape[0])
        for itr in indexes_tr:
            i = itr[0]
            j = itr[1]
            l = itr[2]
            k1 = kmin + (i+1) * info_d['deltak']
            k2 = kmin + (j+1) * info_d['deltak']
            k3 = kmin + (l+1) * info_d['deltak']

            "counting triplets"
            counts[i,j,l] = np.einsum('i,i,i',Jkx[i],Jkx[j],Jkx[l]) / info_d['Ngrid'] ** 3.

            "if there are triangles.. measure them maderfac!! XD"
            if counts[i,j,l] > 0.:

                triangles.append([k1,k2,k3])
                triangles_eff.append(info_d['keff'][[i,j,l]])

                bisp_ijl = np.einsum('i,i,i',Ikx[i],Ikx[j],Ikx[l]) / info_d['Ngrid'] ** 3./info_d['I33']

                p0k_i.append(p0k_msn[i])
                p0k_j.append(p0k_msn[j])
                p0k_l.append(p0k_msn[l])

                b123_arr.append(bisp_ijl/counts[i,j,l])
                bksn_arr.append((p0k_msn[i] + p0k_msn[j] + p0k_msn[l])/ng + 1./ng**2.)
                b123_msn_arr.append(b123_arr[-1] - bksn_arr[-1])

                bkijl[i,j,l] = b123_msn_arr[-1]

                q123_arr.append(bkijl[i,j,l] /(p0k_msn[i]*p0k_msn[j] +
                                               p0k_msn[j]*p0k_msn[l] + p0k_msn[l]*p0k_msn[i]))

    except:

        for i in range(Nmaxbk):
            if i%3 ==0: print("i %d"%i)
            k1 = kmin + (i+1) * info_d['deltak']
            for j in range(i,Nmaxbk):
                k2 = kmin + (j+1) * info_d['deltak']
                for l in range(j,Nmaxbk):
                    k3 = kmin + (l+1) * info_d['deltak']
                    if (k3<(k1+k2)) and ((info_d['deltak'] *  3. + k1 + k2 + k3) < (3. * info_d['kny_pkbktk'][1])): # factor of 3 to account for each k oscillating up to k_i + dk/2

                        "counting triplets"
                        counts[i,j,l] = np.einsum('i,i,i',Jkx[i],Jkx[j],Jkx[l]) / info_d['Ngrid'] ** 3.

                        "if there are triangles.. measure them maderfac!! XD"
                        if counts[i,j,l] > 0.:

                            triangles.append([k1,k2,k3])
                            indexes_tr.append([i,j,l])
                            triangles_eff.append(info_d['keff'][[i,j,l]])

                            bisp_ijl = np.einsum('i,i,i',Ikx[i],Ikx[j],Ikx[l]) / info_d['Ngrid'] ** 3./info_d['I33']

                            p0k_i.append(p0k_msn[i])
                            p0k_j.append(p0k_msn[j])
                            p0k_l.append(p0k_msn[l])

                            b123_arr.append(bisp_ijl/counts[i,j,l])
                            bksn_arr.append((p0k_msn[i] + p0k_msn[j] + p0k_msn[l])/ng + 1./ng**2.)
                            b123_msn_arr.append(b123_arr[-1] - bksn_arr[-1])

                            bkijl[i,j,l] = b123_msn_arr[-1]

                            q123_arr.append(bkijl[i,j,l] /(p0k_msn[i]*p0k_msn[j] +
                                                           p0k_msn[j]*p0k_msn[l] + p0k_msn[l]*p0k_msn[i]))
        np.savez(tr_config_file_name,indexes_tr=np.array(indexes_tr))

    triangles = np.array(triangles)
    print("found and measured %d triangles!"%triangles.shape[0])

    out_bk = {}

    out_bk['xbk']       = np.linspace(1,triangles.shape[0],triangles.shape[0])
    out_bk['triangles'] = triangles
    out_bk['tr_co_eff'] = np.array(triangles_eff)
    out_bk['idx_tr']    = np.array(indexes_tr)
    out_bk['pk_tr']     = np.vstack((np.array(p0k_i),np.array(p0k_j),np.array(p0k_l))).T
    out_bk['b123']      = np.array(b123_arr)
    out_bk['bksn']      = np.array(bksn_arr)
    out_bk['b123_msn']  = np.array(b123_msn_arr)
    out_bk['bkijl']     = bkijl
    out_bk['q123']      = np.array(q123_arr)
    out_bk['counts_bk'] = np.array(counts)

    return out_bk

