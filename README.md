# medir_el_medible
measure power spectrum, bispectrum from simulated boxes

This code at the moment is released “as is”,  I plan to revise it and make it officially public and user-friendly, but I am not there yet.  Call it alpha version.  You are welcome to use it, however at this point be aware that this is not the final product. I recommend you do not yet treat it as a black box for any official (i.e. to be published) material.

If you use this code be sure to read and acknowledge the following papers:

https://arxiv.org/abs/1407.5668
https://arxiv.org/abs/2009.02290


##################################

The code requires FFTW http://fftw.org/

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

In case you have any questions, feel free to drop me an email at dgualdi@icc.ub.edu 


