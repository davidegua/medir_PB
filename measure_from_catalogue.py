import numpy as np
from estimators.pkbktk_estimators import Pk_periodic, Bk_periodic
from estimators.load_input import save_measu, load_data_and_cstep_manyrel
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Ngrid", help="ngrid used to assign the particles to the box",
                    type=int)
parser.add_argument("snap", help="snapshot of the simulations, basically the redshift",
                    type=int)
parser.add_argument("baseindex", help="starting point of the counter",
                    type=int)
parser.add_argument("irel", help="which realisation of the simulation",
                    type=int)
parser.add_argument("dkfac", help="multiplying factor of the standard delta_k",
                    type=float)
parser.add_argument("Nlacing", help="number of interlacing steps",
                    type=int)

parser.add_argument("meas_bk", help="'True' # measure bk",
                    type=int)
parser.add_argument("frac", help="random selected fraction of particles to use",
                    type=float)
parser.add_argument("which_cosmo", help="which set of sims with non-fiducial cosmology", type=str)

args = parser.parse_args()

Ngrid    = args.Ngrid
dkfac    = args.dkfac
snap     = args.snap
irel     = args.irel
Nlacing  = args.Nlacing
meas_bk  = np.bool(args.meas_bk)
frac     = args.frac
which_co = args.which_cosmo

sim_idx = args.baseindex + irel
if which_co=='fid':
    sdetails = 'Quijote_' + which_co + '_'
    with open('local_path.txt', 'r') as file:
        local_path = file.read().replace('\n', '')
    path = local_path + "%d/snapdir_00%d/"%(sim_idx,snap)

else:
    sdetails = 'Quijote_' + which_co + '_'
    local_path = "/projects/QUIJOTE/Snapshots/" + which_co + "/"
    path = local_path + "%d/snapdir_00%d/" % (sim_idx, snap)


print("loading the delta field")
cstep_d, Ik_x, Jk_x = load_data_and_cstep_manyrel(Ngrid,dkfac,snap,path,
                                                  Nlacing=Nlacing,sdetails=sdetails, fraction=frac)
"measure pk"
out_pk = Pk_periodic(cstep_d,Ik_x,Jk_x)

"measure bk"
out_bk = 0
if meas_bk == True: out_bk = Bk_periodic(cstep_d,0.0,cstep_d['kmaxbk'],out_pk['p0k_msn'],Ik_x,Jk_x)

"save data for comparison with analytical model"
save_measu(out_pk,out_bk,Ngrid,cstep_d['Lbox'],cstep_d['s_deta'],sim_idx,snap,dkfac,
           Nlacing,meas_bk,frac)

print("end of the code, measurements stored for realisation # %d"%sim_idx)