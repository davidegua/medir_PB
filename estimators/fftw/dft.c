#include "dft.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <complex.h>
#include <fftw3.h>

// function taken from Hector Gil-Marin's Rustico code and adapted to communicate with Python

extern void r2c_dx_to_dk(double in[], double deltak_re[], double deltak_im[], int ngrid, double L1, double L2, int mode_mass_ass)
{
    fftw_complex *out;
    fftw_plan p;
    int i,j,k;
    double cx,cy,cz;
    double *kx;
    double Pi=(4.*atan(1.));
    long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)

    long int index2;
    kx=malloc(sizeof(double)*(ngrid));
    for(i=0;i<ngrid;i++)
    {
        if(i<ngrid/2+1)
        {kx[i]=i*1.0*(2.0*Pi/(L2-L1));}
        else
        {kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));}
    }

    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(ngridtotr2c));
    fftw_plan_with_nthreads(omp_get_max_threads());
    p =  fftw_plan_dft_r2c_3d(ngrid,ngrid,ngrid,in,out,FFTW_ESTIMATE);

    fftw_execute(p);//FFT
    fftw_destroy_plan(p);

    #pragma omp parallel for private(i,j,k,cx,cy,cz,index2) shared(kx,deltak_re,deltak_im,out,ngrid,ngridtotr2c,mode_mass_ass,Pi)
    for(index2=0;index2<ngridtotr2c;index2++)
    {

        i=(int)(index2/(ngrid*ngrid/2+ngrid));
        j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
        k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

        cx=sin( kx[i]*Pi/(2.*kx[ngrid/2]) )/( kx[i]*Pi/(2.*kx[ngrid/2]) );
        cy=sin( kx[j]*Pi/(2.*kx[ngrid/2]) )/( kx[j]*Pi/(2.*kx[ngrid/2]) );
        cz=sin( kx[k]*Pi/(2.*kx[ngrid/2]) )/( kx[k]*Pi/(2.*kx[ngrid/2]) );
        if(kx[i]==0 || mode_mass_ass==0){cx=1.;}
        if(kx[j]==0 || mode_mass_ass==0){cy=1.;}
        if(kx[k]==0 || mode_mass_ass==0){cz=1.;}

        deltak_re[index2]=creal(out[index2])*pow(cx*cy*cz,- mode_mass_ass*1.);
        deltak_im[index2]=cimag(out[index2])*pow(cx*cy*cz,- mode_mass_ass*1.);
    }

    fftw_free(out);
    free(kx);
    fftw_cleanup();
//in not modified
//deltak_re filled
//deltak_im filled
//out is internally created and destroyed
}
//

// version using the shell assignation algorithm used in Rustico
extern void c2r_dk_to_Ikx_v2(double *deltak_re, double *deltak_im, double **Ikout, long int *lar, long int *IJK,// long int *idL,
                             int ngrid, int nbmax, double Deltak, double kf, double *keffar, int *Neffar)
{
    long int jdx;
    fftw_complex *indki;
    double *outIki;
    int i,j,k,ijk;
    long int ngridr2c = (pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)
    long int ngridtot = pow(ngrid,3);
    long int  line_test, N_eff1, K_eff1_now,K_eff1_bef, N_eff1_unique, l_new;
    double kbin, K_eff1;

    indki  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(ngridr2c));
    outIki = (double*) malloc(sizeof(double)*ngridtot);

    fftw_plan_with_nthreads(omp_get_max_threads());
    fftw_plan p;

    for (int b=0; b<nbmax; b++)
    {
        printf("start loop b %d ",b);
        memset(indki,  0, ngridr2c*sizeof(fftw_complex));//set to 0 elements of in_bis
        memset(outIki, 0, ngridtot*sizeof(double));//set to 0 elements of in_bis

        kbin = Deltak * (b + 1);
//        kbin = Deltak * (b);
        ijk  = pow(kbin/kf,2);
        line_test=(long int)(4.4132*pow(ijk,1.4956));

        if( IJK[line_test]>=ijk)
            {do{line_test=line_test-1;}while(IJK[line_test]>=ijk);}
        else
            {line_test=line_test-1;
             do{line_test=line_test+1;}while(IJK[line_test+1]<ijk);}

        line_test=line_test+1, N_eff1=0, K_eff1=0, N_eff1_unique=0;
//        printf("check point 1,  ");

        do{
            i=(int)(lar[line_test]/(ngrid*ngrid*1.));
            j=(int)((lar[line_test]-i*(ngrid*ngrid))/(ngrid*1.));
            k=lar[line_test]-i*ngrid*ngrid-j*ngrid;
            if(k<=ngrid/2)
            {
                l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;
                __real__ indki[l_new]=deltak_re[l_new];
                __imag__ indki[l_new]=deltak_im[l_new];

                K_eff1+=sqrt(IJK[line_test]*1.)*kf;

                if(N_eff1==0)//First iteration
                    {K_eff1_now=IJK[line_test];
                     K_eff1_bef=-1;}
                else
                    {K_eff1_bef=K_eff1_now;
                     K_eff1_now=IJK[line_test];
                     if(K_eff1_bef!=K_eff1_now){N_eff1_unique++;}}
                N_eff1++;
            }
            else
            {
                k=ngrid-k;
                i=ngrid-i;
                j=ngrid-j;
                if(i==ngrid){i=0;}
                if(j==ngrid){j=0;}
                if(k==ngrid){k=0;}
                l_new=(pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j;
                __real__ indki[l_new]=deltak_re[l_new];//creal
                __imag__ indki[l_new]=deltak_im[l_new];//imag
                K_eff1+=sqrt(IJK[line_test]*1.)*kf;

                if(N_eff1==0)//First iteration
                {   K_eff1_now=IJK[line_test];
                    K_eff1_bef=-1;
                }
                else
                {   K_eff1_bef=K_eff1_now;
                    K_eff1_now=IJK[line_test];

                    if(K_eff1_bef!=K_eff1_now){N_eff1_unique++;}
                }
                 N_eff1++;
            }
            line_test++;
//        }while( sqrt(IJK[line_test])*kf<(0.5 * Deltak + kbin ));
        }while( sqrt(IJK[line_test])*kf<( Deltak + kbin ));

        K_eff1=K_eff1/(N_eff1*1.);

        p =  fftw_plan_dft_c2r_3d(ngrid,ngrid,ngrid,indki,outIki,FFTW_ESTIMATE);

        fftw_execute(p);//FFT

        // store the output of the shell inverse Fourier transform
        for (jdx=0; jdx < ngridtot; jdx++)
            {Ikout[b][jdx] = outIki[jdx];}

        keffar[b] = K_eff1;
        Neffar[b] = N_eff1;

        printf(" --> end loop b %d \n",b);
    }

    fftw_destroy_plan(p);
    fftw_free(indki);
    free(outIki);
    fftw_cleanup();
}

// compute the delta Fourier transform weighted by the legendre polynomial (LOS = z-axis)
extern void r2c_dx_to_dk_multip(double in[], double deltak_req[], double deltak_imq[],double deltak_reh[], double deltak_imh[],
                                int ngrid, double L1, double L2, int mode_mass_ass)
{
    fftw_complex *out;
    fftw_plan p;
    int i,j,k;
    double cx,cy,cz;
    double *kx;
    double Pi=(4.*atan(1.));
    long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)
    double musq, Leg2, Leg4;

    long int index2;
    kx=malloc(sizeof(double)*(ngrid));
    for(i=0;i<ngrid;i++)
    {
        if(i<ngrid/2+1)
        {kx[i]=i*1.0*(2.0*Pi/(L2-L1));}
        else
        {kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));}
    }

    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(ngridtotr2c));
    fftw_plan_with_nthreads(omp_get_max_threads());
    p =  fftw_plan_dft_r2c_3d(ngrid,ngrid,ngrid,in,out,FFTW_ESTIMATE);

    fftw_execute(p);//FFT
    fftw_destroy_plan(p);

    #pragma omp parallel for private(i,j,k,cx,cy,cz,index2,musq,Leg2,Leg4) shared(kx,deltak_req,deltak_imq,deltak_reh,deltak_imh,out,ngrid,ngridtotr2c,mode_mass_ass,Pi)
    for(index2=0;index2<ngridtotr2c;index2++)
    {

        i=(int)(index2/(ngrid*ngrid/2+ngrid));
        j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
        k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

        cx=sin( kx[i]*Pi/(2.*kx[ngrid/2]) )/( kx[i]*Pi/(2.*kx[ngrid/2]) );
        cy=sin( kx[j]*Pi/(2.*kx[ngrid/2]) )/( kx[j]*Pi/(2.*kx[ngrid/2]) );
        cz=sin( kx[k]*Pi/(2.*kx[ngrid/2]) )/( kx[k]*Pi/(2.*kx[ngrid/2]) );
        if(kx[i]==0 || mode_mass_ass==0){cx=1.;}
        if(kx[j]==0 || mode_mass_ass==0){cy=1.;}
        if(kx[k]==0 || mode_mass_ass==0){cz=1.;}

        musq = pow(kx[k],2.)/(pow(kx[i],2.) + pow(kx[j],2.) + pow(kx[k],2.));
        Leg2 = 0.5   * (3.  * musq - 1.);
        Leg4 = 0.125 * (35. * pow(musq,2.) - 30. * musq + 3.);

        deltak_req[index2]=creal(out[index2])*pow(cx*cy*cz,- mode_mass_ass*1.) * Leg2;
        deltak_imq[index2]=cimag(out[index2])*pow(cx*cy*cz,- mode_mass_ass*1.) * Leg2;

        deltak_reh[index2]=creal(out[index2])*pow(cx*cy*cz,- mode_mass_ass*1.) * Leg4;
        deltak_imh[index2]=cimag(out[index2])*pow(cx*cy*cz,- mode_mass_ass*1.) * Leg4;
    }

    fftw_free(out);
    free(kx);
    fftw_cleanup();
//in not modified
//deltak_re filled
//deltak_im filled
//out is internally created and destroyed
}
//
