#include <stdio.h>

/*
Written by JRM, changed for use with GA and simplex codes by RP
*/
void trid(double *,double*, double*,double*, int);

/*
 Given model parameters and upstream boundary conditions, solve transport eqn with transient storage
// cbound is numtim length vector of upstream boundary conditions
// cdown is the downstream values after approximation
// completed is always 0.0 as long as the function returns properly
// *unused is for compatibility with discus
-RP
*/

void approximate_(double *cbound, double *U,double *D,int *numbox, int *numtimptr, double *dxptr, double *dtptr, double  *k1ptr, double *k3ptr,double * unused,double * completed, double *cdown){
double a[1000];
double b[1000];
double c[1000];
double d[1000];
double concentration[1000];
double deadzone[1000];
int n = *numbox;
int numtim = *numtimptr;
double dx = *dxptr, dt = *dtptr, k1 = *k1ptr, k3 = *k3ptr;
int i = 0;
double upstream[5000];
double downstream[5000];
  double theta = 0.5;
// changing this so that it reflects literature and so it matches DISCUS
// Russ had k2 = k1*A/As in CN and k2 = As with A = 1.0 in DISCUS
// this way, input k3 is As like in DISCUS, A=1.0, and k2 can equal k1*A/As in discretization below -RP
  double k2 = 0;
	if(k3 > 0.0) k2 = k1/k3;
  double velocity = *U;
  double dispersion = *D;
  double Courant;
  double Diffusion;
  double alpha;
  double beta;
  double gamma;
  double time = 0.0;
  
// Discretization variables
  alpha = 1.0 + dt*theta*k2;
  beta = dt*k2*theta;
  Courant = velocity*dt/dx;
  Diffusion = dispersion*dt/(dx*dx);
   for(i=0;i<n+1;++i)
   {
    concentration[i] = 0.0;
    deadzone[i] = 0.0;
   }
// because JRM uses 1 based indexing, copy cbound[] to upstream[] 
upstream[0] = 0.0;
   for(i=1;i<numtim+1;++i) // RP
   {
    upstream[i] = cbound[i-1]; // RP
   }
// loop through time
int j =0;
   while(j<numtim) // RP
   {
    // JRM time = time + dt;
	j++; // RP; note that input cbound must be 0 according to JRM
    a[0] = 0.0;
    b[0] = 0.0; 
    c[0] = 0.0;
    d[0] =  0.0;
    a[1] = 0.0;
    b[1] = 1.0; 
    c[1] = 0.0;
    d[1] =  upstream[j];// JRM (int)time];
// Setup LHS
   for(i=2;i<n;++i)
   {
    gamma = deadzone[i]+dt*k2*(1.0-theta)*(concentration[i]-deadzone[i]);
    a[i] = -0.5*theta*Courant - theta*Diffusion;
    b[i] = 1.0 + 2.0*theta*Diffusion + dt*theta*k1*(1.0-beta/alpha); 
    c[i] = 0.5*theta*Courant - theta*Diffusion;
    d[i] =  concentration[i] + k1*dt*theta*gamma/alpha
            - Courant*(1.0-theta)*0.5*(concentration[i+1]-concentration[i-1]) 
            + Diffusion*(1.0-theta)*(concentration[i+1]-2.0*concentration[i]+concentration[i-1])
            - k1*dt*(1.0-theta)*(concentration[i]-deadzone[i]);
   }

    a[n] = -1.0;
    b[n] = 1.0; 
    c[n] = 0.0;
    d[n] =  0.0;

    trid(a,b,c,d,n);

// Do transient storage
   for(i=1;i<n+1;++i)
   {
    gamma = deadzone[i]+dt*k2*(1.0-theta)*(concentration[i]-deadzone[i]);
    concentration[i] = d[i];
    deadzone[i] = (beta/alpha)*concentration[i]+(gamma/alpha);
   }
// copy 
    downstream[j] = concentration[n];
   }
// RP All done time integration

// copy downstream 
   for(i=0;i<numtim;++i)
   {
	if(i<numtim)
		cdown[i] = downstream[i+1]; //RP   
    }

// RP maintain compatibility with discus_only_v2
*completed = 0.0;
return;
} // end

// JRM trid
void trid(double sub[],double diag[],double sup[],double rhs[],int n)
      {
      int i;
/*    
      The tridiagonal linear system
 
      sub[i] * x[i-1] + diag[i] * x[i] + sup[i] * x[i+1] = b[i]  i=1,n
 
      with (sub[1] taken to be zero and sup[n] taken to be zero)
      is solved by factorization and substitution.  

      The factorization is returned in sub, diag, sup and the solution 
      is returned in b.
*/

      if(n<2)
        {
        rhs[1] = rhs[1] / diag[1];
        return;
        }

      for(i=2;i<n+1;++i)
        {
        sub[i] = sub[i] / diag[i-1] ;
        diag[i] = diag[i] - sub[i]*sup[i-1] ;
        rhs[i] = rhs[i] - sub[i]*rhs[i-1] ;
        }
       
        rhs[n] = rhs[n] / diag[n] ;

      for(i=n-1;i>0;--i)
        {
        rhs[i] = (rhs[i] - sup[i]*rhs[i+1]) / diag[i] ;
        }

        return;
 
      }

