#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

complex get_staples(int x, int y, int z, int t, int dir);
void update();
double plaq();
double plaq2();
void printfC(char *s, complex z);
float magz(complex z);
void load_staples(int dir);
int read_ascii_lat(char *fname);
void fluxplane(int z_xyplane);
void update_flux(int z_xyplane);
void meas_flux();
void dumplat();
void dumpplaq();

#define COLD 0
#define HOT 1
#define LOADFILE 2


typedef struct {
   int x;
   int y;
   int z;
   int t;
} parity;


// Global Variables
///////////////////
complex *****L;
complex *****S;
double beta;
int Nvol;
int Np; // number of EVEN, ODD sites: 1/2 Nvol
int Nx; 
int Ny;
int Nz;
int Ns; // N_space
int Nt; //N_time
parity *EVEN;
parity *ODD;
double **mflux;


///////////////////////////////////
// main()
///////////////////////////////////
int main(int argc, char **argv) {
   int x,y,z,t;
   int dir, dir1, dir2;
   int i,j,k,m, warms, trajecs, measure;
   double plaquette;
   double c,s;
   complex stp;
   char loadfile[100] = "";
   int args=argc;
   complex staple;
   //   int Npar;
   int INIT=0;
   int siteno;
   long int seed;
   double r;


   /////////////////////////////////
   // Default parameter values
   /////////////////////////////////
   Ns = 8;
   Nx = 8;
   Ny = 8;
   Nz = 8;
   Nt = 8;
   beta = 1.0;
   warms = 0;
   trajecs = 1;
   measure = 1;
   seed = time(0);


   ////////////////////////////////
   // Parse command line
   ////////////////////////////////

   for(;args>0; --args) {
      if(!(strcmp(argv[args-1], "-u") && strcmp(argv[args-1], "-usage"))) {
         printf("%s command line options (defaults):\n", argv[0]);
	 printf("\t-N (%d) -Ns (%d) -Nx/Ny/Nz -Nt (%d)\n", Ns, Nx, Nt);
         printf("\t-beta (%f)\n", beta);
         printf("\t-warms (%d)\n", warms);
	 printf("\t-trajecs (%d)\n", trajecs);
         printf("\t-meas (%d)\n", measure);
         printf("\t-init [cold, hot, or file filename] (cold)\n");
	 printf("\t-seed (%ld)\n", seed);
         exit(0);
      }
      if(!(strcmp(argv[args-1], "N=") && strcmp(argv[args-1], "-N")))
         Nx = Ny = Nz = Nt = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "Ns=") && strcmp(argv[args-1], "-Ns")))
         Nx = Ny = Nz = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "Nx=") && strcmp(argv[args-1], "-Nx")))
         Nx = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "Ny=") && strcmp(argv[args-1], "-Ny")))
         Ny = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "Nz=") && strcmp(argv[args-1], "-Nz")))
         Nz = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "Nt=") && strcmp(argv[args-1], "-Nt")))
         Nt = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "beta=") && strcmp(argv[args-1], "-beta")))
         beta = atof(argv[args]);
      if(!(strcmp(argv[args-1], "warms=") && strcmp(argv[args-1], "-warms")))
         warms = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "trajecs=") && strcmp(argv[args-1], "-trajecs")))
         trajecs = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "meas=") && strcmp(argv[args-1], "-meas")))
         measure = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "seed=") && strcmp(argv[args-1], "-seed")))
         seed = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "init=") && strcmp(argv[args-1], "-init"))) {
	 printf("init [%s]\n", argv[args]);
	 printf("init [%s] [%s]\n", argv[args], argv[args+1]);
	 if(!strcmp(argv[args], "cold")) { 
	       INIT = COLD;
	 } else if(!strcmp(argv[args], "hot")) { 
	       INIT = HOT;
	 } else if(!strcmp(argv[args], "file")) { 
	       INIT = LOADFILE;
	       strcpy(loadfile, argv[args+1]);
	 }
      }
   }

   Nvol = Nx*Ny*Nz*Nt;
   
   printf("# Nx Ny Nz Nt Nvol = %d %d %d %d %d\n", Nx,Ny,Nz,Nt,Nvol);
   printf("# beta    = %f\n", beta);
   printf("# warms   = %d\n", warms);
   printf("# trajecs = %d\n", trajecs);
   printf("# measure = %d\n", measure);
   printf("# seed    = %d\n", seed);
   if(INIT == COLD) printf("# init = %d cold\n", INIT);
   if(INIT == HOT) printf("# init = %d hot\n", INIT);
   if(INIT == LOADFILE) printf("# init = %d file: %s\n", INIT, loadfile);


   /////////////////////////////////////
   // Initialize
   /////////////////////////////////////

   // seed drand48
   srand48(seed);




   L = (complex *****)malloc(Nt*sizeof(complex ****));
   for(i=0; i<Nt; i++) {
      L[i] = (complex ****)malloc(Nz*sizeof(complex ***));
      for(j=0; j<Nz; j++) {
	 L[i][j] = (complex ***)malloc(Ny*sizeof(complex **));
	 for(k=0; k<Ny; k++) {
	    L[i][j][k] = (complex **)malloc(Nx*sizeof(complex *));
	    for(m=0; m<Nx; m++) {
	       L[i][j][k][m] = (complex *)malloc(4*sizeof(complex));
	    }
	 }
      }
   }

   S = (complex *****)malloc(Nt*sizeof(complex ****));
   for(i=0; i<Nt; i++) {
      S[i] = (complex ****)malloc(Nz*sizeof(complex ***));
      for(j=0; j<Nz; j++) {
	 S[i][j] = (complex ***)malloc(Ny*sizeof(complex **));
	 for(k=0; k<Ny; k++) {
	    S[i][j][k] = (complex **)malloc(Nx*sizeof(complex *));
	    for(m=0; m<Nx; m++) {
	       S[i][j][k][m] = (complex *)malloc(4*sizeof(complex));
	    }
	 }
      }
   }


   // make parity structures
   Np = Nvol/2; 
   EVEN = (parity *)malloc(Np * sizeof(parity));
   ODD = (parity *)malloc(Np * sizeof(parity));

   i=j=0;
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
	       //	       printf("P %d %d %d %d : %i : ",x,y,z,t,i);
	       if(((x+y+z+t)%2)==0) {
		  EVEN[i].x = x;
		  EVEN[i].y = y;
		  EVEN[i].z = z;
		  EVEN[i].t = t;
		  // printf("EVEN\n");
		  i++;
	       } else {
		  ODD[j].x = x;
		  ODD[j].y = y;
		  ODD[j].z = z;
		  ODD[j].t = t;
		  //printf("ODD\n");
		  j++;
	       }
	    }}}}

   
   /*
     for(j=0; j<Np; j++) {
	 x = ODD[j].x;
	 y = ODD[j].y;
	 z = ODD[j].z;
	 t = ODD[j].t; 
	 printf("ODD[ %d ] = (%d %d %d %d)\n",j,x,y,z,t);
     }
     for(i=0; i<Np; i++) {
	 x = EVEN[i].x;
	 y = EVEN[i].y;
	 z = EVEN[i].z;
	 t = EVEN[i].t; 
	 printf("EVEN[ %d ] = (%d %d %d %d)\n",i,x,y,z,t);
     }
   */


   ////////////////////////////////////////
   // make flux measurement array

   mflux = (double **)malloc(Nt*sizeof(double *));
   for(i=0; i<Nt; i++) {
      mflux[i] = (double *)malloc(Nz*sizeof(double));
   }
			   



   ///////////////////////////////////////////////////////////
   // Initialize Lattice
   ///////////////////////////////////////////////////////////


   // cold lattice
   if(INIT == COLD) {
      for(t=0; t<Nt; t++) {
	 for(z=0; z<Nz; z++) {
	    for(y=0; y<Ny; y++) {
	       for(x=0; x<Nx; x++) {
		  for(dir=0; dir<4; dir++) {
		     L[t][z][y][x][dir] = 1.0 + 0.0*I;
		  }}}}}
      printf("# initialized COLD lattice\n");
   }

   // hot lattice
   if(INIT == HOT) {
      for(t=0; t<Nt; t++) {
	 for(z=0; z<Nz; z++) {
	    for(y=0; y<Ny; y++) {
	       for(x=0; x<Nx; x++) {
		  for(dir=0; dir<4; dir++) {
		     r = 4.0*M_PI*(drand48() - 0.5);
		     c = cos(r);
		     s = sin(r);
		     L[t][z][y][x][dir] = (c + s*I);

		  }}}}}
      printf("# initialized HOT lattice\n");
   }


   // Read an ASCII lattice in
   if(INIT == LOADFILE) {
      read_ascii_lat(loadfile);
   }

   ////////////////////////////////////////////////
   // Specialty configs for testing
   ////////////////////////////////////////////////

   // Special lattice to test staples
   /**
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  siteno = Nx*(Ny*(Nz*(t) + z) + y) + x;
		  // L[t][z][y][x][dir] = 1.0 + 0.0*I;
		  L[t][z][y][x][dir] = dir+1 + (((float)dir+1.0)/4)*I;
	       }}}}}
   //exit(0);
   **/

   /**
   // EVEN
   for(dir=0; dir<4; dir++) {
      for(t=0; t<Nt; t++) {
	 for(z=0; z<Nz; z++) {
	    for(y=0; y<Ny; y++) {
	       for(x=0; x<Nx; x++) {
		  if(((x+y+z+t)%2)==0) {
		     staple = get_staples(x,y,z,t,dir);
		     printf("x= %d y= %d z= %d t= %d dir= %d ", x,y,z,t,dir);
		     printf("staple= %f %f\n", creal(staple), -cimag(staple));
		  }}}}}
   }
   // ODD
   for(dir=0; dir<4; dir++) {
      for(t=0; t<Nt; t++) {
	 for(z=0; z<Nz; z++) {
	    for(y=0; y<Ny; y++) {
	       for(x=0; x<Nx; x++) {
		  if(((x+y+z+t)%2)==1) {
		     staple = get_staples(x,y,z,t,dir);
		     printf("x= %d y= %d z= %d t= %d dir= %d ", x,y,z,t,dir);
		     printf("staple= %f %f\n", creal(staple), -cimag(staple));
		  }}}}}
   }
   **/


   ///////////////////////////
   // add a flux plane
   ///////////////////////////
   fluxplane(0);
   printf("# added xy flux plane at z=0\n");
   //   dumplat();
   //   meas_flux(Nz, Nt);
   //   dumpplaq();

   //////////////////////////////////////////////////////
   plaquette = plaq2();
   printf("# initial plaqette %f %f\n", beta, plaquette);
   //////////////////////////////////////////////////////

   //// END PREPARATION ////



   //////////////////////////////////////////////////////
   // do equilibration sweeps
   //////////////////////////////////////////////////////

   for(i=0; i<warms; i++) {
      update_flux(0);
   }
   printf("# %d warms completed\n", warms);
   printf("# starting measurements\n");
   printf("# ------------------------------------\n");

   ///////////////////////////////////////////////////
   // do measurements
   ///////////////////////////////////////////////////
   for(i=0; i<trajecs; i++) {
      update_flux(0);
      if((i % measure) == 0) {
	 plaquette = plaq2();
	 printf("PLAQ %f %f\n", beta, plaquette);
	 meas_flux(Nz, Nt);
      }
   }

   //   printf("FINAL\n");
   //   dumplat();


   return 0;
}




/////////////////////////////////////////////////////
// Functions
/////////////////////////////////////////////////////

void printfC(char *s, complex z) {
   printf("%s %f %f\n", s, creal(z), cimag(z));
} 

float magz(complex z) {
   return sqrt(creal((conj(z)*z)));
}


void dumplat() {
   int x,y,z,t,dir;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
	       printf("%d %d %d %d\n",x,y,z,t); 
	       for(dir=0; dir<4; dir++) {
		  printf("   %d : %f %f t: %f\n", dir,
			 creal(L[t][z][y][x][dir]),
			 cimag(L[t][z][y][x][dir]),
			 atan2(cimag(L[t][z][y][x][dir]),creal(L[t][z][y][x][dir])));

	       }
	    }
	 }
      }   
   }
}

void dumpplaq() {
   complex staple, stptmp, P, p;
   int x,y,z,t,dir,dir1;
   int count;
   complex P01,P02,P03,P12,P13,P23;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
      // p01, p02, p03
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][0]);
      stptmp *= conj(L[t][z][y][x][1]);
      p += stptmp;
      P01 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][0]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      P02 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][0]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P03 = stptmp;
      //    count++;

      // p12, p13
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][1]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      P12 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][1]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P13 = stptmp;
      //    count++;

      // p23
      stptmp  = L[t][z][y][x][2];
      stptmp *= L[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][2]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P23 = stptmp;
      //    count++;

      printf("%d %d %d %d P01 %f %f %f P02 %f %f P03 %f %f P12 %f %f P13 %f %f P23 %f %f\n",
	     x,y,z,t, 
	     creal(P01),cimag(P01), atan2(cimag(P01),creal(P01)),
	     creal(P02),cimag(P02),
	     creal(P03),cimag(P03), creal(P12),cimag(P12),
	     creal(P13),cimag(P13), creal(P23),cimag(P23));
	     

      /*
      printf("%d %d %d %d : (%f,%f) (%f,%f) (%f,%f)\n\t\t(%f,%f) (%f,%f) (%f,%f)\n",
	     x,y,z,t, 
	     creal(P01),cimag(P01), creal(P02),cimag(P02),
	     creal(P03),cimag(P03), creal(P12),cimag(P12),
	     creal(P13),cimag(P13), creal(P23),cimag(P23));
      */
   }}}}

   //   P = p;
   //   p /= Nx*Nx*Nx*Nt*6;
   //   printf("P %f %f, p %f %f\n", creal(P), cimag(P), creal(p), cimag(p));  
   //   printf("count %d\n", count);
   //   return creal(p);
}




void update() {
   int x,y,z,t,dir, i;
   double r, th, scale, oldaction, newaction, actiondiff;
   complex staples, U, Unew, dU;
   double latmag, accrej, pold, pnew, pdiff, accinc, accyes;

   scale = 1.4;

   latmag = 0.0;
   accrej = 0.0;
   //   printf("update: beta= %f\n", beta);

   for(dir=0; dir<4; dir++) {
      load_staples(dir);

      // ODD sites
      for(i=0; i<Np; i++) {

	 x = ODD[i].x;
	 y = ODD[i].y;
	 z = ODD[i].z;
	 t = ODD[i].t; 
	 
	 U = L[t][z][y][x][dir];

	 //	 printf("update %d %d %d %d -> %d p 1 drand: %f\n", 
	 //		x,y,z,t,dir, drand48());

	 th = scale * 2 *(drand48() - 0.5);
	 //	 do{ th = scale * (drand48() - 0.5); }
	 //	 while( fabs(th) < 0.25*scale );
	 //	 c = cos(theta); s = sin(theta);

	 //	 printf("th= %f\n", th);
	 dU = cos(th) + I*sin(th);
	 Unew = dU * U;
	 Unew /= magz(Unew);  // reunitarize
	 //	 staples = get_staples(x,y,z,t,dir);
	 staples = S[t][z][y][x][dir];

	 oldaction = beta * creal(U * staples);
	 newaction = beta * creal(Unew * staples);
	 
	 //  actiondiff = newaction - oldaction;
	 
	 if(newaction > oldaction) {
	    L[t][z][y][x][dir] = Unew;
	    accrej++; accyes=1;
	    //printf("accept: new > old\n");
	 } else {
	    r = drand48();
	    if( r < exp(newaction - oldaction)) {
	       L[t][z][y][x][dir] = Unew;
	       accrej++; //accyes=1;
	       //printf("accept: exp(new-old) > r %f\n", r);
	    }
	 }
      }
   }

   for(dir=0; dir<4; dir++) {
      load_staples(dir);

      // EVEN sites
      for(i=0; i<Np; i++) {
	 x = EVEN[i].x;
	 y = EVEN[i].y;
	 z = EVEN[i].z;
	 t = EVEN[i].t;
	 
	 U = L[t][z][y][x][dir];
	 //	 printf("update %d %d %d %d -> %d p 2 drand: %f\n", 
	 //		x,y,z,t,dir,drand48());
	 th = scale * 2 * (drand48() - 0.5);
	 //	 printf("th= %f\n", th);
	 dU = cos(th) + I*sin(th);
	 Unew = dU * U;
	 Unew /= magz(Unew);  // reunitarize
	 // staples = get_staples(x,y,z,t,dir);
	 staples = S[t][z][y][x][dir];
	 
	 oldaction = beta * creal(U * staples);
	 newaction = beta * creal(Unew * staples);
	 
	 //  actiondiff = newaction - oldaction;
	 
	 if(newaction > oldaction) {
	    L[t][z][y][x][dir] = Unew;
	    accrej++; accyes=1;
	    //printf("accept: new > old\n");
	 } else {
	    r = drand48();
	    if( r < exp(newaction - oldaction)) {
	       L[t][z][y][x][dir] = Unew;
	       accrej++; //accyes=1;
	       //printf("accept: exp(new-old) > r %f\n", r);
	    }
	 }
      }
   }
   //   accrej /= Nvol*4;
   //   printf("# acc/rej= %f\n", accrej);
}
		    

/////////////////////////////////////////
// Puts flux w through all XY planes
/////////////////////////////////////////

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


//////////////////////////////////////////
// Do NOT update [0] links on z_xyplane 
//////////////////////////////////////////

void update_flux(int z_xyplane) {
   int x,y,z,t,dir, i;
   double r, th, scale, oldaction, newaction, actiondiff;
   complex staples, U, Unew, dU;
   double latmag, accrej, pold, pnew, pdiff, accinc, accyes;

   scale = 1.4;

   latmag = 0.0;
   accrej = 0.0;


   for(dir=0; dir<4; dir++) {
      load_staples(dir);

      // ODD sites
      for(i=0; i<Np; i++) {

	 x = ODD[i].x;
	 y = ODD[i].y;
	 z = ODD[i].z;
	 t = ODD[i].t; 

	 if((z==z_xyplane) && ((dir==0) || (dir==1))) continue;

	 
	 U = L[t][z][y][x][dir];

	 //	 printf("update %d %d %d %d -> %d p 1 drand: %f\n", 
	 //		x,y,z,t,dir, drand48());

	 th = scale * 2 *(drand48() - 0.5);
	 //	 do{ th = scale * (drand48() - 0.5); }
	 //	 while( fabs(th) < 0.25*scale );
	 //	 c = cos(theta); s = sin(theta);

	 //	 printf("th= %f\n", th);
	 dU = cos(th) + I*sin(th);
	 Unew = dU * U;
	 Unew /= magz(Unew);  // reunitarize
	 //	 staples = get_staples(x,y,z,t,dir);
	 staples = S[t][z][y][x][dir];

	 oldaction = beta * creal(U * staples);
	 newaction = beta * creal(Unew * staples);
	 
	 //  actiondiff = newaction - oldaction;
	 
	 if(newaction > oldaction) {
	    L[t][z][y][x][dir] = Unew;
	    accrej++; accyes=1;
	    //printf("accept: new > old\n");
	 } else {
	    r = drand48();
	    if( r < exp(newaction - oldaction)) {
	       L[t][z][y][x][dir] = Unew;
	       accrej++; //accyes=1;
	       //printf("accept: exp(new-old) > r %f\n", r);
	    }
	 }
      }
   }

   for(dir=0; dir<4; dir++) {
      load_staples(dir);

      // EVEN sites
      for(i=0; i<Np; i++) {
	 x = EVEN[i].x;
	 y = EVEN[i].y;
	 z = EVEN[i].z;
	 t = EVEN[i].t;

	 if((z==z_xyplane) && ((dir==0) || (dir==1))) continue;

	 
	 U = L[t][z][y][x][dir];
	 //	 printf("update %d %d %d %d -> %d p 2 drand: %f\n", 
	 //		x,y,z,t,dir,drand48());
	 th = scale * 2 * (drand48() - 0.5);
	 //	 printf("th= %f\n", th);
	 dU = cos(th) + I*sin(th);
	 Unew = dU * U;
	 Unew /= magz(Unew);  // reunitarize
	 // staples = get_staples(x,y,z,t,dir);
	 staples = S[t][z][y][x][dir];
	 
	 oldaction = beta * creal(U * staples);
	 newaction = beta * creal(Unew * staples);
	 
	 //  actiondiff = newaction - oldaction;
	 
	 if(newaction > oldaction) {
	    L[t][z][y][x][dir] = Unew;
	    accrej++; accyes=1;
	    //printf("accept: new > old\n");
	 } else {
	    r = drand48();
	    if( r < exp(newaction - oldaction)) {
	       L[t][z][y][x][dir] = Unew;
	       accrej++; //accyes=1;
	       //printf("accept: exp(new-old) > r %f\n", r);
	    }
	 }
      }
   }
   //   accrej /= Nvol*4;
   //   printf("# acc/rej= %f\n", accrej);
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





double plaq() {
   int x,y,z,t,dir;
   double plaq;
   complex p, U, staple;

   plaq = 0.0;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  staple = get_staples(x,y,z,t,dir);
		  U = L[t][z][y][x][dir];
		  //p = U * conj(staple);
		  p = U * staple;
		  plaq += creal(p);
   }}}}}
   plaq /= Nvol*24;

   return plaq;
}

void load_staples(int dir) {
   int x,y,z,t;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
		  S[t][z][y][x][dir] = get_staples(x,y,z,t,dir);
   }}}}
}


complex get_staples(int x, int y, int z, int t, int dir) {
   complex staple, stptmp;
   
   switch(dir) {
   case 0:
      // upward staples
      stptmp  = L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple = stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][z][(y-1+Nx)%Nx][(x+1)%Nx][1]);
      stptmp *= conj(L[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= L[t][z][(y-1+Nx)%Nx][x][1];
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nx)%Nx][y][(x+1)%Nx][2]);
      stptmp *= conj(L[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= L[t][(z-1+Nx)%Nx][y][x][2];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][y][(x+1)%Nx][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 1:
      // upward staples
      stptmp  = L[t][z][(y+1)%Nx][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = L[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][z][(y+1)%Nx][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nx)%Nx][(y+1)%Nx][x][2]);
      stptmp *= conj(L[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= L[t][(z-1+Nx)%Nx][y][x][2];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][(y+1)%Nx][x][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 2:
      // upward staples
      stptmp  = L[t][(z+1)%Nx][y][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[t][(z+1)%Nx][y][x][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = L[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][(z+1)%Nx][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[t][(z+1)%Nx][(y-1+Nx)%Nx][x][1]);
      stptmp *= conj(L[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= L[t][z][(y-1+Nx)%Nx][x][1];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][(z+1)%Nx][y][x][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 3:
      // upward staples
      stptmp  = L[(t+1)%Nt][z][y][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[(t+1)%Nt][z][y][x][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = L[(t+1)%Nt][z][y][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[(t+1)%Nt][z][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][z][(y-1+Nx)%Nx][x][1]);
      stptmp *= conj(L[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= L[t][z][(y-1+Nx)%Nx][x][1];
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][(z-1+Nx)%Nx][y][x][2]);
      stptmp *= conj(L[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= L[t][(z-1+Nx)%Nx][y][x][2];
      staple += stptmp;
      break;
   }
    
   return staple;
}

double plaq2() {
   complex staple, stptmp, P, p;
   int x,y,z,t,dir,dir1;
   int count;
   complex P01,P02,P03,P12,P13,P23;

   p = 0.0 + 0.0*I;
   // count = 0;
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {

      // p01, p02, p03
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][0]);
      stptmp *= conj(L[t][z][y][x][1]);
      p += stptmp;
      P01 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][0]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      P02 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][0]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P03 = stptmp;
      //    count++;

      // p12, p13
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][1]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      P12 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][1]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P13 = stptmp;
      //    count++;

      // p23
      stptmp  = L[t][z][y][x][2];
      stptmp *= L[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][2]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P23 = stptmp;
      //    count++;

      /*
      printf("%d %d %d %d : (%f,%f) (%f,%f) (%f,%f)\n\t\t(%f,%f) (%f,%f) (%f,%f)\n",
	     x,y,z,t, 
	     creal(P01),cimag(P01), creal(P02),cimag(P02),
	     creal(P03),cimag(P03), creal(P12),cimag(P12),
	     creal(P13),cimag(P13), creal(P23),cimag(P23));
      */
   }}}}

   p /= Nvol*6;
   //   printf("P %f %f, p %f %f\n", creal(P), cimag(P), creal(p), cimag(p));  
   //   printf("count %d\n", count);
   return creal(p);
}



int read_ascii_lat(char *fname) {
   char line[100];
   FILE *fp;
   int x,y,z,t,i,a,dir;
   float Ur, Ui;

   i = 0;
   fp = fopen(fname, "r");

   printf("reading from %s\n", fname);

   // discard header
   while (i<3) {
      if( fgets (line, 100 , fp)!=NULL ) {
	 i++;
      }
   }

   for(t=0; t<Nt; t++)
      for(z=0; z<Nz; z++)
	 for(y=0; y<Ny; y++)
	    for(x=0; x<Nx; x++)
	       for(dir=0; dir<4; dir++) {
		  // get x= y= z= t= line
		     if (fgets(line, 150, fp) == NULL) break;
                  //read a link
		  for(i=0; i<9; i++) {
		     fscanf(fp,"%f %f\n", &Ur, &Ui);
		     if(i==0) {
			printf("coords %d %d %d %d dir %d Ur: %f Ui: %f\n", x,y,z,t,dir,Ur, Ui);
			L[t][z][y][x][dir] = Ur + I*Ui;
		     }
		  }
	       }

   /*
   for(dir=0; dir<4; dir++) {
      printf("%f %f\n", 
	     creal(L[11][11][11][11][dir]), 
	     cimag(L[11][11][11][11][dir]));
   }
   */

   fclose(fp);

}
		  
