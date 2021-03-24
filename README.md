# medir_el_medible
measure power spectrum, bispectrum and trispectrum from simulated boxes


I acknowledge the 
pySpectrum (https://github.com/changhoonhahn/pySpectrum)
Rustico (https://github.com/hectorgil/Rustico/)
Pylians3 (https://github.com/franciscovillaescusa/Pylians3) 

codes from which i have used many useful functions (cited in the code).

Run with

Arguments:
 Ngrid    = 128
 snap     = 1 
 base_rel = 8000 #starting index for the realisations
 irel     = 1
 dkfac    = 3    #dk/kf 
 Nlacing  = 2    # N interlacing steps
 tksym    = all  # 'all' for ... configurations, 'asym' for ... configurations (no unconnected term)

python measure_from_catalogue.py  128 1 8000 1 3 1 all
