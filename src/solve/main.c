#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
/* 
by Rick Page
A wrapper for either DISCUS or central difference approximation (OTIS/CN)

1D solver for advection diffusion transient storage. 
Given U,D,k1,k2, as well as domain information M,N,dx,dt, this program reads M upstream concentrations from a 
space-delimited file named lb.dat and records downstream concentrations into rb.dat. 
Usage: 

Input:
M - number of observations
N - number of grid points in space
U - advective velocity
D - diffusion coefficient
k1- transient storage coefficient
k2- ratio of main channel to transient storage cross sectional area
dt- timestep
dx- spacing along open-channel

Output:
The downstream values are written to the file rb.dat
 
Example:
This example simulates a stream with upstream concentration 1.0 at each of ten intervals. The length 
of the reach is 2.5 meters, so N is set to 25 and dx to 0.1. There is transient storage so k1 and k2 must be
set to non-zero values, in this case 0.01 and 0.1. Velocity of the flow is 0.5 m/s and the dispesion
coefficient is set to 0.05 m^2/s:
> echo "1 1 1 1 1 1 1 1 1 1" > lb.dat
> ./bin/DISCUS 10 25 0.5 0.05 0.01 0.1 1 .1
> cat rb.dat 
0.000000 0.000004 0.000665 0.027044 0.245139 0.613679 0.844978 0.935266 0.968190 0.981448

Note:
in this file, approximate_ is the name of the C and the FORTRAN routines for either method, testtridJRM.c and discus_only_v2.f respectively, so only one C code is needed to call either routine. 

*/

// Declare external solver; depending on which .o has been linked, this will be CN or DISCUS
extern void approximate_(double *cbound, double *u,double *d,int *numbox, int *numtim,double *xdelta, double* tdelta, double * k1, double *k2,double *d1,double * dk2, double *cdown);

// Helper functions
int readData(double *C, int N, char *filename){
	int i;
	FILE * f;
		
	if(!filename) { 
		fprintf(stderr,"Error: filename not provided"); 
		return 0; }
	f = fopen(filename,"r");
	if(!f) { 
		fprintf(stderr,"Error: %s cannot be opened for read.",filename); 
		return 0; }
	if(!C) { fprintf(stderr,"Error: *C points to nothing: allocate memory first.(%s)",__func__);
		return 0;		
		}
	for(i = 0; i < N; i++){
		if(fscanf(f,"%lf",&C[i]) == EOF) {
			fprintf(stderr,"Error: Could only read %d of %d total mesh points"
			" from file %s(user supplied N too large?)",i,N,filename);
		exit(2);
		}
	}
	fclose(f);
	return 1;
}


int writeData(double *C, int N, char *filename){
	int i;
	FILE * f;
	if(!filename) { 
		fprintf(stderr,"Error: filename not provided"); 
		return 0; }
	f = fopen(filename,"w");
	if(!f) { 
		fprintf(stderr,"Error: %s cannot be opened for write.",filename); 
		return 0; }
	for(i = 0; i < N; i++){
		fprintf(f,"%f ",C[i]);	
	}
	fclose(f);
	return 1;
}

/*
Memory
*/
// Allocate memory, and display informative message on failure 
double * safeMalloc(size_t size,int line, char *file){
	double * ptr = malloc(size);
	if(!ptr){
		fprintf(stderr,"Error: Failed to allocate %lu bytes. line: %d , %s",size,line,file);
		exit(1);
	}
	return ptr;
}
//
// MAIN
//
int main(int argc, char **argv){
	double d1 = 0.0, dk2 = 0.0;
	double dx=1.0,dt=1.0,k2=0.0,k1=0.0,D=0.0,U=1.0;
int T,N;

// INPUT
	if(argc < 9){
		printf("Solver for 1D advection diffusion transient storage.\n\
Given U,D,k1,k2, as well as domain information M,N,dx,dt, this program reads M upstream concentrations from a \n\
space-delimited file named lb.dat and records downstream concentrations into rb.dat. \n\n\
Usage:\n%s M N U D k1 k2 dt dx\nreads lb.dat as upstream and records rb.dat as downstream", argv[0]);
printf("1D solver for advection diffusion transient storage. \n\
Given U,D,k1,k2, as well as domain information M,N,dx,dt, this program reads M upstream\n\
concentrations from a space-delimited file named lb.dat and records downstream \n\
concentrations into rb.dat.\n\n\
Usage: \n\
%s M N U D k1 k2 dt dx\n\n\
Input:\n\
M - number of observations\n\
N - number of grid points in space\n\
U - advective velocity\n\
D - diffusion coefficient\n\
k1- transient storage coefficient\n\
k2- ratio of main channel to transient storage cross sectional area\n\
dt- timestep\n\
dx- spacing along open-channel\n\n\
Output:\n\
The downstream values are written to the file rb.dat.\n\n\
Example:\n\
This example simulates a stream with upstream concentration 1.0 at each of ten intervals.  \n\
The length of the reach is 2.5 meters, so N is set to 25 and dx to 0.1. There is transient \n\
storage so k1 and k2 must be set to non-zero values, in this case 0.01 and 0.1. Velocity  \n\
of the flow is 0.5 m/s and the dispesion coefficient is set to 0.05 m^2/s: \n\
> echo \"1 1 1 1 1 1 1 1 1 1\" > lb.dat \n\
> ./bin/DISCUS 10 25 0.5 0.05 0.01 0.1 1 .1 \n\
> cat rb.dat  \n\
0.000000 0.000004 0.000665 0.027044 0.245139 0.613679 0.844978 0.935266 0.968190 0.981448 \n\n",argv[0]);

		exit(EXIT_FAILURE);
	} else {
		switch(argc){
			default:
			case 9:
			dx = atof(argv[8]);
			case 8:
			dt = atof(argv[7]);
			case 7:
			k2 = atof(argv[6]);
			case 6:
			k1 = atof(argv[5]);
			case 5:
			D = atof(argv[4]);
			case 4:
			T = atoi(argv[1]);
			N = atoi(argv[2]);

			if(N > INT_MAX) {
				fprintf(stderr,"Input Error: N(%d) > INT_MAX(%d) is not allowed.",N,INT_MAX);
				exit(1);
			} 
			U = atof(argv[3]);
			break;
		}
	}

	size_t	size = (T)*sizeof(double);
	int steps = T-1;
	
	if(steps <= 0){ 
		fprintf(stderr,"\nError: < 1 steps being used."); 
		exit(1);
	}
	// Allocate memory for upstream and downstream
	double * cbound = safeMalloc(size,__LINE__,__FILE__),
		* cdown = safeMalloc(size,__LINE__,__FILE__);
	
	// populate upstream vector with information from lb.dat
	readData(cbound,T,"lb.dat");
	// run solver
	approximate_((double *)cbound, &U,&D,&N,&steps,&dx,&dt,&k1,&k2,&d1,&dk2,(double *)cdown);
	// dk2 must be 0, or the solver has failed
	if(dk2 > 0.0){ 
		fprintf(stderr,"rb.dat not written because computation became unstable"); 	
	}
	else {
		// If the solver succeeded, print rb.dat as the downstream concentrations
		writeData(cdown,T,"rb.dat");	
	}	
	// Free allocated memory
	free(cbound);
	free(cdown);

return dk2;
}
