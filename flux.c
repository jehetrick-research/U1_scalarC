/////////////////////////////////////////
// U(1) Magnetic Fields
/////////////////////////////////////////


//////////////////////////////////////////////////
// Puts z-directed flux w through all XY planes
//////////////////////////////////////////////////

#include "u1includes.h"

// v.1.5
// Assumes COLD lattice first
void fluxplane_twist(int z_xyplane) {
   int x,y,z,t,dir,q=1;


   // Put the Flux plane just one z=0 
   z = z_xyplane;
   for(t=0; t<Nt; t++) {
      //   for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
	       L[t][z][y][x][0] = cos(w*y) - I*sin(w*y);
	    }
	 }
	 //   }
   }   
}	       

void meas_flux_twist() {
   complex stptmp, P;
   int x,y,z,t,dir,dir1;
   int i,j, q;
   complex P01; //,P02,P03,P12,P13,P23;
   complex sumReP01, wphase;
   double th, sumth, thexptn;

   for(i=0; i<Nt; i++)
      for(j=0; j<Nz; j++) 
	 mflux[i][j] = 0.0;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 sumReP01 = 0 + I*0;
	 sumth = 0;
	 for(y=0; y<Ny; y++) { 
	    wphase = (y==Ny-1) ? cos(Ny*w) - I*sin(Ny*w) : 1.0 + I*0.0;
	    for(x=0; x<Nx; x++) {
	       // P01
	       stptmp  = L[t][z][y][x][0];
	       stptmp *= L[t][z][y][(x+1)%Nx][1];
	       stptmp *= conj(wphase * L[t][z][(y+1)%Ny][x][0]);
	       stptmp *= conj(L[t][z][y][x][1]);
	       P01 = stptmp;
	       th = atan2(cimag(P01),creal(P01));
	       /**
		  printf("%d %d %d %d : (%f,%f) th= %f\n",
		  x,y,z,t, creal(P01),cimag(P01), th); 
	       **/
	       
	       sumReP01 += P01;
	       sumth += th;
	    }
	 }
	 mflux[t][z] = sumth/Nx/Ny;
      }}

   for(j=0; j<Nz; j++) {
      for(i=1; i<Nt; i++) {
	 mflux[0][j] += mflux[i][j];
      }
   }

   
   printf("FLUX %f ", beta);
   thexptn = 0.0;
   for(j=0; j<Nz; j++) {
      mflux[0][j] /= Nt;
      //      printf("%f ", atan2( creal(mflux[0][j]), mflux[0][j]) );
      printf("%f ", mflux[0][j] );
      if(j!=0) thexptn += mflux[0][j];
   }
   if(w != 0) printf(": %f %f", thexptn/(Nz-1)/w, mflux[0][Nz/2]/w); else printf("w=0");
   printf("\n");

}





// v.1

void fluxplane(int z_xyplane) {
   int x,y,z,t,dir,q=1;
   double w;

   q = 1;
   w = 2*M_PI*q/(Nx);

   //   z = 0;
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
	       L[t][z][y][x][0] = cos(w*y) - I*sin(w*y);
	       //	       L[t][z][y][x][1] = cos(w*x) + I*sin(w*x);
	    }
	 }
      }
   }   
}	       

void meas_flux() {
   complex stptmp, P;
   int x,y,z,t,dir,dir1;
   int i,j, q;
   complex P01; //,P02,P03,P12,P13,P23;
   complex sumReP01, Pthref;
   double th, sumth, thexptn;

   q = 1;
   Pthref = 2*M_PI*q/Nx;
   
   for(i=0; i<Nt; i++)
      for(j=0; j<Nz; j++) 
	 mflux[i][j] = 0.0;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 sumReP01 = 0 + I*0;
	 sumth = 0;

	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
      // p01, p02, p03
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][0]);
      stptmp *= conj(L[t][z][y][x][1]);
      P01 = stptmp;
      th = atan2(cimag(P01),creal(P01));
      /**
      printf("%d %d %d %d : (%f,%f) th= %f\n",
	     x,y,z,t, creal(P01),cimag(P01), th); 
      **/

      sumReP01 += P01;
      sumth += th;
	    }}
	 //	 mflux[t][z] = sumReP01/Nx/Nx;
	 mflux[t][z] = sumth/Nx/Nx;
      }}

   for(j=0; j<Nz; j++) {
      for(i=1; i<Nt; i++) {
	 mflux[0][j] += mflux[i][j];
      }
   }

   
   printf("FLUX %f ", beta);
   thexptn = 0.0;
   for(j=0; j<Nz; j++) {
      mflux[0][j] /= Nt;
      //      printf("%f ", atan2( creal(mflux[0][j]), mflux[0][j]) );
      printf("%f ", mflux[0][j] );
      if(j!=0) thexptn += mflux[0][j];
   }
   printf(": %f", thexptn/(Nz-1)/Pthref);
   printf("\n");

}


