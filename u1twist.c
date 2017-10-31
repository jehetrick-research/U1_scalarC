////////////////////////////////////////////////////
// 4-dimensional U(1) lattice gauge theory
//
// J.Hetrick Oct.--- 2016
//////////////////////////////////////////////////// 



#include "u1includes.h"


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
double w;
int zplane=0;

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
   w = 0.0;

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
	 printf("\t-w (%f)\n", w);
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
      if(!(strcmp(argv[args-1], "w=") && strcmp(argv[args-1], "-w")))
         w = atof(argv[args]);
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
   printf("# w flux  = %f\n", w);
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



   /////////////////////////////////////////////
   // add a flux plane @ z=zplane
   /////////////////////////////////////////////
   fluxplane_twist(zplane);
   printf("# added xy flux plane at z=0\n");
   //dumplat();

   meas_flux_twist(Nz, Nt);
   //dumpplaq_twist();


   /////////////////////////////////////////////////////
   // Initial Plaquette
   /////////////////////////////////////////////////////
   plaquette = plaq_twist();
   printf("# initial plaqette %f %f\n", beta, plaquette);
   //////////////////////////////////////////////////////

   //// END PREPARATION ////



   // Staples testing v. milc code
   ///////////////////////////////
   /**
   for(dir=0; dir<4; dir++) {
     for(j=0; j<Np; j++) {
	 x = EVEN[j].x;
	 y = EVEN[j].y;
	 z = EVEN[j].z;
	 t = EVEN[j].t; 
	 staple = get_staples_twist(x,y,z,t,dir);
	 printf("L %d %d %d %d [%d]: %f %f\n",x,y,z,t,dir,
		creal(L[t][z][y][x][dir]), cimag(L[t][z][y][x][dir]) );
	 printf("S %d %d %d %d [%d]: %f %f\n",x,y,z,t,dir,
		creal(staple), cimag(staple) );
     }
     for(i=0; i<Np; i++) {
	 x = ODD[i].x;
	 y = ODD[i].y;
	 z = ODD[i].z;
	 t = ODD[i].t; 
	 staple = get_staples_twist(x,y,z,t,dir);
	 printf("L %d %d %d %d [%d]: %f %f\n",x,y,z,t,dir,
		creal(L[t][z][y][x][dir]), cimag(L[t][z][y][x][dir]) );
	 printf("S %d %d %d %d [%d]: %f %f\n",x,y,z,t,dir,
		creal(staple), cimag(staple) );
     }
   }
   **/




   //////////////////////////////////////////////////////
   // do equilibration sweeps
   //////////////////////////////////////////////////////

   for(i=0; i<warms; i++) {
      update_flux(zplane);
   }
   printf("# %d warms completed\n", warms);
   printf("# starting measurements\n");
   printf("# ------------------------------------\n");

   ///////////////////////////////////////////////////
   // do measurements
   ///////////////////////////////////////////////////
   for(i=0; i<trajecs; i++) {
      update_flux(zplane);
      if((i % measure) == 0) {
	 plaquette = plaq_twist();
	 printf("PLAQ %f %f\n", beta, plaquette);
	 meas_flux_twist(Nz, Nt);
      }
   }

   //   printf("FINAL\n");
   //   dumplat();


   return 0;
}
//// END MAIN() /////////////////////////////////////
/////////////////////////////////////////////////////





		    





