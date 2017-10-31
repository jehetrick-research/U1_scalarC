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
complex polyakov_P();


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
complex *****U;
complex *****S;
complex ***P;
double beta;
int Nvol;
int Np;
int Nx;
int Nt;
parity *EVEN;
parity *ODD;


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
   double r;
   long int seed;
   complex Ploop;

   Nx = 8;
   Nt = 8;
   beta = 1.0;
   warms = 0;
   trajecs = 1;
   measure = 1;
   seed = time(0);

   for(;args>0; --args) {
      if(!(strcmp(argv[args-1], "-h") && strcmp(argv[args-1], "--help"))) {
         printf("%s command line options (defaults):\n", argv[0]);
	 printf("\t-N (%d)  -Nx (%d) -Nt (%d)\n", Nx, Nx, Nt);
	 printf("\t or: -L 8x8\n");
         printf("\t-beta (%f)\n", beta);
         printf("\t-warms (%d)\n", warms);
	 printf("\t-trajecs (%d)\n", trajecs);
         printf("\t-meas (%d)\n", measure);
	 printf("\t-seed (%ld)\n", seed);
         printf("\t-init [cold, hot, or file filename] (cold)\n");
         exit(0);
      }
      if(!(strcmp(argv[args-1], "N=") && strcmp(argv[args-1], "-N")))
         Nx = Nt = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "Nx=") && strcmp(argv[args-1], "-Nx")))
         Nx = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "Nt=") && strcmp(argv[args-1], "-Nt")))
         Nt = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "-L")))
	 sscanf(argv[args],"%ix%i", &Nx, &Nt);
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

   printf("# Nx Nt = %d %d\n", Nx, Nt);
   printf("# beta    = %f\n", beta);
   printf("# warms   = %d\n", warms);
   printf("# trajecs = %d\n", trajecs);
   printf("# measure = %d\n", measure);
   printf("# seed    = %ld\n", seed);
   if(INIT == COLD) printf("# init = %d cold\n", INIT);
   if(INIT == HOT) printf("# init = %d hot\n", INIT);
   if(INIT == LOADFILE) printf("# init = %d file: %s\n", INIT, loadfile);
   printf("# MEAS_FORMAT: beta plaq re(P) im(P) |P|\n");


   // seed drand48 with the time 
   // srand48(time(0));
   srand48(seed);


   // create Lattice: U_mu(x) array
   U = (complex *****)malloc(Nt*sizeof(complex ****));
   for(i=0; i<Nt; i++) {
      U[i] = (complex ****)malloc(Nx*sizeof(complex ***));
      for(j=0; j<Nx; j++) {
	 U[i][j] = (complex ***)malloc(Nx*sizeof(complex **));
	 for(k=0; k<Nx; k++) {
	    U[i][j][k] = (complex **)malloc(Nx*sizeof(complex *));
	    for(m=0; m<Nx; m++) {
	       U[i][j][k][m] = (complex *)malloc(4*sizeof(complex));
	    }
	 }
      }
   }

   // make staple sum array
   S = (complex *****)malloc(Nt*sizeof(complex ****));
   for(i=0; i<Nt; i++) {
      S[i] = (complex ****)malloc(Nx*sizeof(complex ***));
      for(j=0; j<Nx; j++) {
	 S[i][j] = (complex ***)malloc(Nx*sizeof(complex **));
	 for(k=0; k<Nx; k++) {
	    S[i][j][k] = (complex **)malloc(Nx*sizeof(complex *));
	    for(m=0; m<Nx; m++) {
	       S[i][j][k][m] = (complex *)malloc(4*sizeof(complex));
	    }
	 }
      }
   }

// make Polyakov loop array (only need one complex for each spatial site)
   P = (complex ***)malloc(Nx*sizeof(complex **));
   for(i=0; i<Nx; i++) {
      P[i] = (complex **)malloc(Nx*sizeof(complex *));
      for(j=0; j<Nx; j++) {
	 P[i][j] = (complex *)malloc(Nx*sizeof(complex));
      }
   }


   // make parity structures
   Nvol = Nx*Nx*Nx*Nt;
   Np = Nvol/2; 
   EVEN = (parity *)malloc(Np * sizeof(parity));
   ODD = (parity *)malloc(Np * sizeof(parity));

   i=j=0;
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
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





   ///////////////////////////////////////////////////////////
   // Initialize Lattice
   // 
   // Should be done with argv ("cold" "hot" or "filename")
   ///////////////////////////////////////////////////////////


   // cold lattice
   if(INIT == COLD) {
      for(t=0; t<Nt; t++) {
	 for(z=0; z<Nx; z++) {
	    for(y=0; y<Nx; y++) {
	       for(x=0; x<Nx; x++) {
		  for(dir=0; dir<4; dir++) {
		     U[t][z][y][x][dir] = 1.0 + 0.0*I;
		  }}}}}
      printf("# initialized COLD lattice\n");
   }

   // hot lattice
   if(INIT == HOT) {
      for(t=0; t<Nt; t++) {
	 for(z=0; z<Nx; z++) {
	    for(y=0; y<Nx; y++) {
	       for(x=0; x<Nx; x++) {
		  for(dir=0; dir<4; dir++) {
		     c = 2.0*(drand48() - 0.5);
		     s = 2.0*(drand48() - 0.5);
		     U[t][z][y][x][dir] = (c + s*I)/sqrt(c*c + s*s);

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
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  siteno = Nx*(Nx*(Nx*(t) + z) + y) + x;
		  // U[t][z][y][x][dir] = 1.0 + 0.0*I;
		  U[t][z][y][x][dir] = dir+1 + (((float)dir+1.0)/4)*I;
	       }}}}}
   //exit(0);
   **/

   /**
   // EVEN
   for(dir=0; dir<4; dir++) {
      for(t=0; t<Nt; t++) {
	 for(z=0; z<Nx; z++) {
	    for(y=0; y<Nx; y++) {
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
	 for(z=0; z<Nx; z++) {
	    for(y=0; y<Nx; y++) {
	       for(x=0; x<Nx; x++) {
		  if(((x+y+z+t)%2)==1) {
		     staple = get_staples(x,y,z,t,dir);
		     printf("x= %d y= %d z= %d t= %d dir= %d ", x,y,z,t,dir);
		     printf("staple= %f %f\n", creal(staple), -cimag(staple));
		  }}}}}
   }
   **/

   //////////////////////////////////////////////////////

   //   plaquette = plaq();
   //   printf("# initial plaqette1 %f %f\n", beta, plaquette);

   plaquette = plaq2();
   Ploop = polyakov_P();
   printf("INIT %f %f %f %f %f\n", beta, plaquette, creal(Ploop), cimag(Ploop), magz(Ploop));
   



   //////////////////////////////////////////////////////
   // do equilibration sweeps
   //////////////////////////////////////////////////////

   for(i=0; i<warms; i++) {
      update();
   }
   printf("# %d warms completed\n", warms);
   printf("# starting measurements\n");
   printf("# ------------------------------------\n");

   ///////////////////////////////////////////////////
   // do measurements
   ///////////////////////////////////////////////////
   for(i=0; i<trajecs; i++) {
      update();
      if((i % measure) == 0) {
	 plaquette = plaq2();
	 Ploop = polyakov_P();
	 printf("MEAS %f %f %f %f %f\n", beta, plaquette, creal(Ploop), cimag(Ploop), magz(Ploop));
      }
   }

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


void update() {
   int x,y,z,t,dir, i;
   double r, th, scale, oldaction, newaction, actiondiff;
   complex staples, Uold, Unew, dU;
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
	 
	 Uold = U[t][z][y][x][dir];

	 //	 printf("update %d %d %d %d -> %d p 1 drand: %f\n", 
	 //		x,y,z,t,dir, drand48());

	 th = scale * 2 *(drand48() - 0.5);
	 //	 do{ th = scale * (drand48() - 0.5); }
	 //	 while( fabs(th) < 0.25*scale );
	 //	 c = cos(theta); s = sin(theta);

	 //	 printf("th= %f\n", th);
	 dU = cos(th) + I*sin(th);
	 Unew = dU * Uold;
	 // Unew /= magz(Unew);  // reunitarize
	 //	 staples = get_staples(x,y,z,t,dir);
	 staples = S[t][z][y][x][dir];

	 oldaction = beta * creal(Uold * staples);
	 newaction = beta * creal(Unew * staples);
	 
	 //  actiondiff = newaction - oldaction;
	 
	 if(newaction > oldaction) {
	    U[t][z][y][x][dir] = Unew;
	    accrej++; accyes=1;
	    //printf("accept: new > old\n");
	 } else {
	    r = drand48();
	    if( r < exp(newaction - oldaction)) {
	       U[t][z][y][x][dir] = Unew;
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
	 
	 Uold = U[t][z][y][x][dir];
	 //	 printf("update %d %d %d %d -> %d p 2 drand: %f\n", 
	 //		x,y,z,t,dir,drand48());
	 th = scale * 2 * (drand48() - 0.5);
	 //	 printf("th= %f\n", th);
	 dU = cos(th) + I*sin(th);
	 Unew = dU * Uold;
	 // Unew /= magz(Unew);  // reunitarize
	 // staples = get_staples(x,y,z,t,dir);
	 staples = S[t][z][y][x][dir];
	 
	 oldaction = beta * creal(Uold * staples);
	 newaction = beta * creal(Unew * staples);
	 
	 //  actiondiff = newaction - oldaction;
	 
	 if(newaction > oldaction) {
	    U[t][z][y][x][dir] = Unew;
	    accrej++; accyes=1;
	    //printf("accept: new > old\n");
	 } else {
	    r = drand48();
	    if( r < exp(newaction - oldaction)) {
	       U[t][z][y][x][dir] = Unew;
	       accrej++; //accyes=1;
	       //printf("accept: exp(new-old) > r %f\n", r);
	    }
	 }
      }
   }
   //   accrej /= Nt*Nx*Nx*Nx*4;
   //   printf("# acc/rej= %f\n", accrej);
}
		    

double plaq() {
   int x,y,z,t,dir;
   double plaq;
   complex p, link, staple;

   plaq = 0.0;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  staple = get_staples(x,y,z,t,dir);
		  link = U[t][z][y][x][dir];
		  //p = link * conj(staple);
		  p = link * staple;
		  plaq += creal(p);
   }}}}}
   plaq /= Nx*Nx*Nx*Nt*24;

   return plaq;
}

void load_staples(int dir) {
   int x,y,z,t;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
		  S[t][z][y][x][dir] = get_staples(x,y,z,t,dir);
   }}}}
}



complex get_staples(int x, int y, int z, int t, int dir) {
   complex staple, stptmp;
   
   switch(dir) {
   case 0:
      // upward staples
      stptmp  = U[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(U[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(U[t][z][y][x][1]);
      staple = stptmp;
      stptmp  = U[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(U[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(U[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = U[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(U[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(U[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(U[t][z][(y-1+Nx)%Nx][(x+1)%Nx][1]);
      stptmp *= conj(U[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= U[t][z][(y-1+Nx)%Nx][x][1];
      staple += stptmp;
      stptmp  = conj(U[t][(z-1+Nx)%Nx][y][(x+1)%Nx][2]);
      stptmp *= conj(U[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= U[t][(z-1+Nx)%Nx][y][x][2];
      staple += stptmp;
      stptmp  = conj(U[(t-1+Nt)%Nt][z][y][(x+1)%Nx][3]);
      stptmp *= conj(U[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= U[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 1:
      // upward staples
      stptmp  = U[t][z][(y+1)%Nx][x][0];
      stptmp *= conj(U[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(U[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = U[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(U[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(U[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = U[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(U[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(U[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(U[t][z][(y+1)%Nx][(x-1+Nx)%Nx][0]);
      stptmp *= conj(U[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= U[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(U[t][(z-1+Nx)%Nx][(y+1)%Nx][x][2]);
      stptmp *= conj(U[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= U[t][(z-1+Nx)%Nx][y][x][2];
      staple += stptmp;
      stptmp  = conj(U[(t-1+Nt)%Nt][z][(y+1)%Nx][x][3]);
      stptmp *= conj(U[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= U[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 2:
      // upward staples
      stptmp  = U[t][(z+1)%Nx][y][x][0];
      stptmp *= conj(U[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(U[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = U[t][(z+1)%Nx][y][x][1];
      stptmp *= conj(U[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(U[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = U[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(U[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(U[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(U[t][(z+1)%Nx][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(U[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= U[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(U[t][(z+1)%Nx][(y-1+Nx)%Nx][x][1]);
      stptmp *= conj(U[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= U[t][z][(y-1+Nx)%Nx][x][1];
      staple += stptmp;
      stptmp  = conj(U[(t-1+Nt)%Nt][(z+1)%Nx][y][x][3]);
      stptmp *= conj(U[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= U[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 3:
      // upward staples
      stptmp  = U[(t+1)%Nt][z][y][x][0];
      stptmp *= conj(U[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(U[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = U[(t+1)%Nt][z][y][x][1];
      stptmp *= conj(U[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(U[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = U[(t+1)%Nt][z][y][x][2];
      stptmp *= conj(U[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(U[t][z][y][x][2]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(U[(t+1)%Nt][z][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(U[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= U[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(U[(t+1)%Nt][z][(y-1+Nx)%Nx][x][1]);
      stptmp *= conj(U[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= U[t][z][(y-1+Nx)%Nx][x][1];
      staple += stptmp;
      stptmp  = conj(U[(t+1)%Nt][(z-1+Nx)%Nx][y][x][2]);
      stptmp *= conj(U[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= U[t][(z-1+Nx)%Nx][y][x][2];
      staple += stptmp;
      break;
   }
    
   return staple;
}

double plaq2() {
   complex staple, stptmp, P, p;
   int x,y,z,t,dir,dir1;
   int count;

   p = 0.0 + 0.0*I;
   // count = 0;
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {

      // p01, p02, p03
      stptmp  = U[t][z][y][x][0];
      stptmp *= U[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(U[t][z][(y+1)%Nx][x][0]);
      stptmp *= conj(U[t][z][y][x][1]);
      p += stptmp;
      //    count++;
      stptmp  = U[t][z][y][x][0];
      stptmp *= U[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(U[t][(z+1)%Nx][y][x][0]);
      stptmp *= conj(U[t][z][y][x][2]);
      p += stptmp;
      //    count++;
      stptmp  = U[t][z][y][x][0];
      stptmp *= U[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(U[(t+1)%Nt][z][y][x][0]);
      stptmp *= conj(U[t][z][y][x][3]);
      p += stptmp;
      //    count++;

      // p12, p13
      stptmp  = U[t][z][y][x][1];
      stptmp *= U[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(U[t][(z+1)%Nx][y][x][1]);
      stptmp *= conj(U[t][z][y][x][2]);
      p += stptmp;
      //    count++;
      stptmp  = U[t][z][y][x][1];
      stptmp *= U[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(U[(t+1)%Nt][z][y][x][1]);
      stptmp *= conj(U[t][z][y][x][3]);
      p += stptmp;
      //    count++;

      // p23
      stptmp  = U[t][z][y][x][2];
      stptmp *= U[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(U[(t+1)%Nt][z][y][x][2]);
      stptmp *= conj(U[t][z][y][x][3]);
      p += stptmp;
      //    count++;
   }

	 }}}
   P = p;
   p /= Nx*Nx*Nx*Nt*6;
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
      for(z=0; z<Nx; z++)
	 for(y=0; y<Nx; y++)
	    for(x=0; x<Nx; x++)
	       for(dir=0; dir<4; dir++) {
		  // get x= y= z= t= line
		  //   if (fgets(line, 150, fp) == NULL) break;
                  //read a link
		  for(i=0; i<9; i++) {
		     fscanf(fp,"%f %f\n", &Ur, &Ui);
		     if(i==0) {
			printf("coords %d %d %d %d dir %d Ur: %f Ui: %f\n", x,y,z,t,dir,Ur, Ui);
			U[t][z][y][x][dir] = Ur + I*Ui;
		     }
		  }
	       }

   /*
   for(dir=0; dir<4; dir++) {
      printf("%f %f\n", 
	     creal(U[11][11][11][11][dir]), 
	     cimag(U[11][11][11][11][dir]));
   }
   */

   fclose(fp);

}
		  

complex polyakov_P() {
   int x,y,z,t;
   complex p, Pave;

   Pave = 0.0;
   for(x=0; x<Nx; x++)
   for(y=0; y<Nx; y++)
   for(z=0; z<Nx; z++) {
      p = U[0][z][y][x][3];
      for(t=1; t<Nt; t++) {
	 p *= U[t][z][y][x][3];
      }
      P[x][y][z] = p;
      Pave += p;
   }

   Pave /= Nx*Nx*Nx;
   return Pave;
}
	   
