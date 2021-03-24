# medir_el_medible
measure power spectrum, bispectrum and trispectrum from simulated boxes

it requires FFTW http://fftw.org/

I acknowledge the 
pySpectrum (https://github.com/changhoonhahn/pySpectrum)
Rustico (https://github.com/hectorgil/Rustico/)
Pylians3 (https://github.com/franciscovillaescusa/Pylians3) 

codes from which i have used many useful functions (cited in the code).

Run with

Arguments:
 Ngrid    = 128
 snap     = 3    #simulation snapshot 
 base_rel = 8000 #starting index for the realisations
 irel     = 0
 dkfac    = 3    #dk/kf 
 Nlacing  = 2    # N interlacing steps
 measBk   = 1    # measure or not the bispectrum
 frac     = 1.   # use 100% of the particles
 cosmo    = fid  # use simulations for a given cosmology (Quijote simulations)

python measure_from_catalogue.py 256 3 0 0 3 2 1 1 fid
