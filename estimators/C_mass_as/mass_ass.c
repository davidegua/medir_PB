#include "mass_ass.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// functions from Rustico code by Hector Gil-Marin

void pq5s_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);

    int i,j,l;
    int  xindex=(int)(ngrid*s_x/(L2-L1));
    int  yindex=(int)(ngrid*s_y/(L2-L1));
    int  zindex=(int)(ngrid*s_z/(L2-L1));

    int  xindex2;
    int  yindex2;
    int  zindex2;
    long int index;
    double sx;
    double sy;
    double sz;

    double factor_x,factor_y,factor_z;
    factor_x=0;
    factor_y=0;
    factor_z=0;
    int number_of_cells_explored=3;

    for(i=-number_of_cells_explored;i<=number_of_cells_explored;i++)
    {
        xindex2=xindex+i;
        if(xindex2<0){xindex2=xindex2+ngrid;}
        if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
        sx=(-(xindex+i+0.5)+s_x*ngrid*1./(L2-L1));
        //
        factor_x=0;
        if( fabs(sx)<1.)//only contributing to the adjacent grids
        {
                  factor_x=11./20.-0.5*sx*sx+0.25*sx*sx*sx*sx-1./12.*pow(fabs(sx),5);
        }
        if( fabs(sx)>=1. && fabs(sx)<2.)//only contributing to the adjacent grids
        {
                  factor_x=17./40.+5./8.*fabs(sx)-7./4.*sx*sx+5./4.*pow(fabs(sx),3)-3./8.*pow(sx,4)+1./24.*pow(fabs(sx),5);
        }
        if( fabs(sx)>=2. && fabs(sx)<3.)//only contributing to the adjacent grids
        {
                  factor_x=243./120.-81./24.*fabs(sx)+9./4.*pow(sx,2)-3./4.*pow(fabs(sx),3)+1./8.*pow(sx,4)-1./120.*pow(fabs(sx),5);
        }

        //
        for(j=-number_of_cells_explored;j<=number_of_cells_explored;j++)
        {
            yindex2=yindex+j;
            if(yindex2<0){yindex2=yindex2+ngrid;}
            if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
            sy=(-(yindex+j+0.5)+s_y*ngrid*1./(L2-L1));
            factor_y=0;
            if( fabs(sy)<1.)//only contributing to the adjacent grids
            {
                      factor_y=11./20.-0.5*sy*sy+0.25*sy*sy*sy*sy-1./12.*pow(fabs(sy),5);
            }
            if( fabs(sy)>=1. && fabs(sy)<2.)//only contributing to the adjacent grids
            {
                      factor_y=17./40.+5./8.*fabs(sy)-7./4.*sy*sy+5./4.*pow(fabs(sy),3)-3./8.*pow(sy,4)+1./24.*pow(fabs(sy),5);
            }
            if( fabs(sy)>=2. && fabs(sy)<3.)//only contributing to the adjacent grids
            {
                      factor_y=243./120.-81./24.*fabs(sy)+9./4.*pow(sy,2)-3./4.*pow(fabs(sy),3)+1./8.*pow(sy,4)-1./120.*pow(fabs(sy),5);
            }

            //
            for(l=-number_of_cells_explored;l<=number_of_cells_explored;l++)
            {
                zindex2=zindex+l;
                if(zindex2<0){zindex2=zindex2+ngrid;}
                if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
                sz=(-(zindex+l+0.5)+s_z*ngrid*1./(L2-L1));
                factor_z=0;
                if( fabs(sz)<1.)//only contributing to the adjacent grids
                {
                          factor_z=11./20.-0.5*sz*sz+0.25*sz*sz*sz*sz-1./12.*pow(fabs(sz),5);
                }
                if( fabs(sz)>=1. && fabs(sz)<2.)//only contributing to the adjacent grids
                {
                          factor_z=17./40.+5./8.*fabs(sz)-7./4.*sz*sz+5./4.*pow(fabs(sz),3)-3./8.*pow(sz,4)+1./24.*pow(fabs(sz),5);
                }
                if( fabs(sz)>=2. && fabs(sz)<3.)//only contributing to the adjacent grids
                {
                          factor_z=243./120.-81./24.*fabs(sz)+9./4.*pow(sz,2)-3./4.*pow(fabs(sz),3)+1./8.*pow(sz,4)-1./120.*pow(fabs(sz),5);
                }
                index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);
                //					  delta[index]=delta[index]+factor_x*factor_y*factor_z*weight;
                delta[index] = factor_x*factor_y*factor_z*weight;

            }
        }
    }

}

extern void mass_ass_pq5s(double *delta, double *sx, double *sy, double *sz, double *weight_arr, double L2, double L1,
                          int ngrid, int Ndata)
{
    int c;
    #pragma omp parallel for private(c) shared(ngrid,Ndata,delta)
    for(c=0; c<Ndata; c++)
    {
//        pq5s_assingment(delta_data, sx[c], sy[c], sz[c], weight[c], L2, L1, ngrid);
	double s_x=sx[c]+(-L1);
	double s_y=sy[c]+(-L1);
	double s_z=sz[c]+(-L1);
    double weight = weight_arr[c];
    int i,j,l;
    int  xindex=(int)(ngrid*s_x/(L2-L1));
    int  yindex=(int)(ngrid*s_y/(L2-L1));
    int  zindex=(int)(ngrid*s_z/(L2-L1));

    int  xindex2;
    int  yindex2;
    int  zindex2;
    long int index;
    double sx;
    double sy;
    double sz;

    double factor_x,factor_y,factor_z;
    factor_x=0;
    factor_y=0;
    factor_z=0;
    int number_of_cells_explored=3;

    for(i=-number_of_cells_explored;i<=number_of_cells_explored;i++)
    {
        xindex2=xindex+i;
        if(xindex2<0){xindex2=xindex2+ngrid;}
        if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
        sx=(-(xindex+i+0.5)+s_x*ngrid*1./(L2-L1));
        //
        factor_x=0;
        if( fabs(sx)<1.)//only contributing to the adjacent grids
        {
                  factor_x=11./20.-0.5*sx*sx+0.25*sx*sx*sx*sx-1./12.*pow(fabs(sx),5);
        }
        if( fabs(sx)>=1. && fabs(sx)<2.)//only contributing to the adjacent grids
        {
                  factor_x=17./40.+5./8.*fabs(sx)-7./4.*sx*sx+5./4.*pow(fabs(sx),3)-3./8.*pow(sx,4)+1./24.*pow(fabs(sx),5);
        }
        if( fabs(sx)>=2. && fabs(sx)<3.)//only contributing to the adjacent grids
        {
                  factor_x=243./120.-81./24.*fabs(sx)+9./4.*pow(sx,2)-3./4.*pow(fabs(sx),3)+1./8.*pow(sx,4)-1./120.*pow(fabs(sx),5);
        }

        //
        for(j=-number_of_cells_explored;j<=number_of_cells_explored;j++)
        {
            yindex2=yindex+j;
            if(yindex2<0){yindex2=yindex2+ngrid;}
            if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
            sy=(-(yindex+j+0.5)+s_y*ngrid*1./(L2-L1));
            factor_y=0;
            if( fabs(sy)<1.)//only contributing to the adjacent grids
            {
                      factor_y=11./20.-0.5*sy*sy+0.25*sy*sy*sy*sy-1./12.*pow(fabs(sy),5);
            }
            if( fabs(sy)>=1. && fabs(sy)<2.)//only contributing to the adjacent grids
            {
                      factor_y=17./40.+5./8.*fabs(sy)-7./4.*sy*sy+5./4.*pow(fabs(sy),3)-3./8.*pow(sy,4)+1./24.*pow(fabs(sy),5);
            }
            if( fabs(sy)>=2. && fabs(sy)<3.)//only contributing to the adjacent grids
            {
                      factor_y=243./120.-81./24.*fabs(sy)+9./4.*pow(sy,2)-3./4.*pow(fabs(sy),3)+1./8.*pow(sy,4)-1./120.*pow(fabs(sy),5);
            }

            //
            for(l=-number_of_cells_explored;l<=number_of_cells_explored;l++)
            {
                zindex2=zindex+l;
                if(zindex2<0){zindex2=zindex2+ngrid;}
                if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
                sz=(-(zindex+l+0.5)+s_z*ngrid*1./(L2-L1));
                factor_z=0;
                if( fabs(sz)<1.)//only contributing to the adjacent grids
                {
                          factor_z=11./20.-0.5*sz*sz+0.25*sz*sz*sz*sz-1./12.*pow(fabs(sz),5);
                }
                if( fabs(sz)>=1. && fabs(sz)<2.)//only contributing to the adjacent grids
                {
                          factor_z=17./40.+5./8.*fabs(sz)-7./4.*sz*sz+5./4.*pow(fabs(sz),3)-3./8.*pow(sz,4)+1./24.*pow(fabs(sz),5);
                }
                if( fabs(sz)>=2. && fabs(sz)<3.)//only contributing to the adjacent grids
                {
                          factor_z=243./120.-81./24.*fabs(sz)+9./4.*pow(sz,2)-3./4.*pow(fabs(sz),3)+1./8.*pow(sz,4)-1./120.*pow(fabs(sz),5);
                }
                index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);
                delta[index]=delta[index]+factor_x*factor_y*factor_z*weight;

            }
        }
     }
     }
}



void pcs_assingment(double delta[], double s_x, double s_y, double s_z, double weight, double L2, double L1, int ngrid)
{
	s_x=s_x+(-L1);
	s_y=s_y+(-L1);
	s_z=s_z+(-L1);


    int i,j,l;
    int  xindex=(int)(ngrid*s_x/(L2-L1));
    int  yindex=(int)(ngrid*s_y/(L2-L1));
    int  zindex=(int)(ngrid*s_z/(L2-L1));
    long int index;
    int  xindex2;
    int  yindex2;
    int  zindex2;

    double sx;
    double sy;
    double sz;

    double factor_x,factor_y,factor_z;
    factor_x=0;
    factor_y=0;
    factor_z=0;
    int number_of_cells_explored=2;

    for(i=-number_of_cells_explored;i<=number_of_cells_explored;i++)
    {
        xindex2=xindex+i;
        if(xindex2<0){xindex2=xindex2+ngrid;}
        if(xindex2>ngrid-1){xindex2=xindex2-ngrid;}
        sx=(-(xindex+i+0.5)+s_x*ngrid*1./(L2-L1));

        factor_x=0;
        if( fabs(sx)<1)//only contributing to the adjacent grids
        {
                  factor_x=2./3.-sx*sx+1./2.*pow(fabs(sx),3);
        }
        if( fabs(sx)>=1 && fabs(sx)<2.)//only contributing to the adjacent grids
        {
                  factor_x=4./3.-2.*fabs(sx)+sx*sx-1./6.*pow(fabs(sx),3);
        }

        for(j=-number_of_cells_explored;j<=number_of_cells_explored;j++)
        {
            yindex2=yindex+j;
            if(yindex2<0){yindex2=yindex2+ngrid;}
            if(yindex2>ngrid-1){yindex2=yindex2-ngrid;}
            sy=(-(yindex+j+0.5)+s_y*ngrid*1./(L2-L1));
            factor_y=0;
            if( fabs(sy)<1)//only contributing to the adjacent grids
            {
                      factor_y=2./3.-sy*sy+1./2.*pow(fabs(sy),3);
            }
            if( fabs(sy)>=1 && fabs(sy)<2.)//only contributing to the adjacent grids
            {
                      factor_y=4./3.-2.*fabs(sy)+sy*sy-1./6.*pow(fabs(sy),3);
            }


            for(l=-number_of_cells_explored;l<=number_of_cells_explored;l++)
            {
                zindex2=zindex+l;
                if(zindex2<0){zindex2=zindex2+ngrid;}
                if(zindex2>ngrid-1){zindex2=zindex2-ngrid;}
                sz=(-(zindex+l+0.5)+s_z*ngrid*1./(L2-L1));
                factor_z=0;
                if( fabs(sz)<1)//only contributing to the adjacent grids
                {
                        factor_z=2./3.-sz*sz+1./2.*pow(fabs(sz),3);
                }
                if( fabs(sz)>=1 && fabs(sz)<2.)//only contributing to the adjacent grids
                {
                        factor_z=4./3.-2.*fabs(sz)+sz*sz-1./6.*pow(fabs(sz),3);
                }
                index=(pow(ngrid,2)*xindex2+ngrid*yindex2+zindex2);
                delta[index]=delta[index]+factor_x*factor_y*factor_z*weight;

			}
	    }
    }
}