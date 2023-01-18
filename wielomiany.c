#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
	
double fi(double x, int k){
	if(k==0){
		return(1);
	}
	if(k==1){
		return(1-x);
	}
	if(k>=2){
        return ((2*(k-1)+1-x)*fi(x,k-1)-(k-1)*fi(x,k-2))/k;
	}
}

/* Pierwsza pochodna fi */	
double dfi(double x, int k){
	if(k==0){
		return(0);
	}
	if(k==1){
		return(-1);
	}
	if(k>=2){
        return (-fi(x,k-1)+(2*(k-1)+1-x)*dfi(x,k-1)-(k-1)*dfi(x,k-2))/k;
	}
}

/* Druga pochodna fi */
double d2fi(double x, int k){
	if(k==0){
		return(0);
	}
	if(k==1){
		return(0);
	}
	if(k>=2){
        return (-2*dfi(x,k-1)+(2*(k-1)+1-x)*d2fi(x,k-1)-(k-1)*d2fi(x,k-2))/k;
	}
}

/* Trzecia pochodna fi */
double d3fi(double x, int k){
	if(k==0){
		return(0);
	}
	if(k==1){
		return(0);
	}
	if(k>=2){	
        return (-3*d2fi(x,k-1)+(2*(k-1)+1-x)*d3fi(x,k-1)-(k-1)*d3fi(x,k-2))/k;
	}
}

void make_spl(points_t * pts, spline_t * spl){
    int n=pts->n;
    double sum;
    double *x=pts->x;
    double *y=pts->y;
    double a = x[0];
    double b = x[n - 1];
    int i, j, k;
    int nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
    char *nbEnv = getenv("APPROX_BASE_SIZE");
    if (nbEnv!=NULL&&atoi(nbEnv)>0)
        nb=atoi(nbEnv);
    matrix_t *eqs=make_matrix(nb, nb + 1);
    for (i=0;i<nb;i++){
        for(j=0;j<nb;j++){
            sum=0;
            for(k=0;k<n;k++)
                sum+=fi(x[k],i)*fi(x[k],j);
                
            put_entry_matrix(eqs, i, j, sum);
        }
    }
    for(i=0;i<nb;i++){
        sum=0;
        for(k=0;k<n;k++)
            sum+=fi(x[k],i)*y[k];

        put_entry_matrix(eqs,i,nb,sum);
    }
    
	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
	
	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi(xx,k);
				spl->f1[i] += ck * dfi(xx,k);
				spl->f2[i] += ck * d2fi(xx,k);
				spl->f3[i] += ck * d3fi(xx,k);
			}
		}
	}
}
