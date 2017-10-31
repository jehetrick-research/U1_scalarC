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
void load_staples();
int read_ascii_lat(char *fname);
int siteno;

#define Nx 12
#define Nt 12

complex L[Nt][Nx][Nx][Nx][4];
complex S[Nt][Nx][Nx][Nx][4];
double beta;


int main(int argc, char **argv) {
   int x,y,z,t;
   int dir, dir1, dir2;
   int i, warms, trajecs, measure;
   double plaquette;
   double c,s;
   complex stp;
   char loadfile[100] = "";
   int args=argc;
   complex staple;

   /*
   if(argc != 5) {
      printf("Usage: %s beta warms trajecs measure\n", argv[0]);
      exit(0);
   } else {
      beta = atof(argv[1]);
      warms = atoi(argv[2]);
      trajecs = atoi(argv[3]);
      measure = atoi(argv[4]);
   }
   */


   //   loadfile = "";
   for(;args>0; --args) {
      if(!(strcmp(argv[args-1], "-u") && strcmp(argv[args-1], "-usage"))) {
         printf("%s command line options (defaults):\n", argv[0]);
         printf("\t-beta (1.0)\n");
         printf("\t-warms (1000)\n");
	 printf("\t-trajecs (100)\n");
         printf("\t-meas (10)\n");
         printf("\t-reload [filename] ()\n");
         exit(0);
      }
      if(!(strcmp(argv[args-1], "beta=") && strcmp(argv[args-1], "-beta")))
         beta = atof(argv[args]);
      if(!(strcmp(argv[args-1], "warms=") && strcmp(argv[args-1], "-warms")))
         warms = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "trajecs=") && strcmp(argv[args-1], "-trajecs")))
         trajecs = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "meas=") && strcmp(argv[args-1], "-meas")))
         measure = atoi(argv[args]);
      if(!(strcmp(argv[args-1], "reload=") && strcmp(argv[args-1], "-reload")))
         strcpy(loadfile, argv[args]);
   }

   printf("# beta    = %f\n", beta);
   printf("# warms   = %d\n", warms);
   printf("# trajecs = %d\n", trajecs);
   printf("# measure = %d\n", measure);
   printf("# reload = %s\n", loadfile);




   // seed drand48 with the time 
   srand48(time(0));


   ///////////////////////////////////////////////
   // Initialize Lattice
   // 
   // Should be done with argv ("cold" "hot" or "filename")

   // cold lattice
   /**
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  L[t][z][y][x][dir] = 1.0 + 0.0*I;
   }}}}}
   **/

   // Read an ASCII lattice in
   read_ascii_lat(loadfile);


   // Special lattice to test staples
   /**
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  siteno = Nx*(Nx*(Nx*(t) + z) + y) + x;
		  // L[t][z][y][x][dir] = 1.0 + 0.0*I;
		  L[t][z][y][x][dir] = dir+1 + (((float)dir+1.0)/4)*I;
	       }}}}}
   **/
   //exit(0);




   /**/
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

   // ODD
      //   for(dir=0; dir<4; dir++) {
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
   /**/



   plaquette = plaq2();
   printf("PLAQ %f %f\n", beta, plaquette);






   // do equilibration
   for(i; i<warms; i++) {
      update();
   }
   printf("# %d warms completed\n", warms);

   for(i=0; i<trajecs; i++) {
      update();
      if((i % measure) == 0) {
	 plaquette = plaq2();
	 printf("PLAQ %f %f\n", beta, plaquette);
      }
   }

   return 0;
}

void printfC(char *s, complex z) {
   printf("%s %f %f\n", s, creal(z), cimag(z));
} 

float magz(complex z) {
   return sqrt(creal((conj(z)*z)));
}


void update() {
   int x,y,z,t,dir;
   double r, th, scale, oldaction, newaction, actiondiff;
   complex staples, U, Unew, dU;
   double latmag, accrej, pold, pnew, pdiff, accinc, accyes;

   scale = 0.5;
   latmag = 0.0;
   accrej = 0.0;
   //   printf("update: beta= %f\n", beta);

   //   load_staples();

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  U = L[t][z][y][x][dir];
		  th = scale * M_PI * 2*(drand48() - 0.5);
		  //printf("th= %f\n", th);
		  dU = cos(th) + I*sin(th);
		  Unew = dU * U;
		  //staples = S[t][z][y][x][dir];
		  staples = get_staples(x,y,z,t,dir);
		  //printfC("U", U);
		  //printfC("dU", dU);
		  //printf("|dU|= %f\n", magz(dU)); 
		  //printfC("Unew", Unew);
		  //printf("|Unew|= %f\n", magz(Unew)); 
		  //printfC("staples=", staples);
		  //printf("|staples|= %f\n", magz(staples)); 
		  //printfC("U*conj(staples)", U * conj(staples));
		  //printfC("Unew*conj(staples)", Unew * conj(staples));

		  //oldaction = beta * creal(U * conj(staples));
		  //newaction = beta * creal(Unew * conj(staples));

		  oldaction = beta * creal(U * staples);
		  newaction = beta * creal(Unew * staples);
		  //printf("noconj action: old %f new %f\n", oldaction, newaction);
		  actiondiff = newaction - oldaction;
		  //printf("action: old %f new %f diff %f\n", 
		  //	 oldaction, newaction, actiondiff);
		  
		  //pold = plaq2();
		  //printf("total beta*Re(plaq_old): %f\n", beta*pold);
		  
		  //accinc = accrej;
		  //accyes=0;
		  if(newaction > oldaction) {
		     L[t][z][y][x][dir] = Unew;
		     accrej++; accyes=1;
		     //printf("accept: new > old\n");
		  } else {
		     r = drand48();
		     //printf("exp(new-old) %f >? r %f\n",exp(newaction-oldaction),r);
		     if( exp(newaction - oldaction) > r) {
			L[t][z][y][x][dir] = Unew;
			accrej++; //accyes=1;
		        //printf("accept: exp(new-old) > r %f\n", r);
		     }
		  }
		  //printf("X= %d %d %d %d dir %d\n",x,y,z,t,dir);
		  //if(accyes == 1) {
		  //pnew = plaq2();
		  //   [M 88printf("total beta*Re(plaq_new): %f\n", beta*pnew);
		  //   printf("DELTA anew - aold: %f\n", actiondiff);
		  //   printf("DELTA pnew - pold: %f\n", pnew - pold);
		  //   //if(abs(actiondiff - (pnew-pold))) {
		  //	printf("DELTAERROR: %f\n", actiondiff - (pnew-pold));
		  //   //}
		  //}

		  //		  latmag += magz(L[t][z][y][x][dir]);
	       }}}}}
   //accrej /= Nt*Nx*Nx*Nx*4;
   //   printf("acc/rej= %f\n", accrej);
}
		    

double plaq() {
   int x,y,z,t,dir;
   double plaq;
   complex p, U, staple;

   plaq = 0.0;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  staple = get_staples(x,y,z,t,dir);
		  U = L[t][z][y][x][dir];
		  p = U * conj(staple);
		  //p = U * staple;
		  plaq += creal(p);
   }}}}}
   plaq /= Nx*Nx*Nx*Nt*24;

   return plaq;
}

void load_staples() {
   int x,y,z,t,dir;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  S[t][z][y][x][dir] = get_staples(x,y,z,t,dir);
   }}}}}
}



complex get_staples(int x, int y, int z, int t, int dir) {
   complex staple, stptmp;
   
   switch(dir) {
   case 0:
      // upward staples
      stptmp  = L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      printf("SF1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple = stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      printf("SF2: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      printf("SF3: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][z][(y-1+Nx)%Nx][(x+1)%Nx][1]);
      stptmp *= conj(L[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= L[t][z][(y-1+Nx)%Nx][x][1];
      printf("SB1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nx)%Nx][y][(x+1)%Nx][2]);
      stptmp *= conj(L[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= L[t][(z-1+Nx)%Nx][y][x][2];
      printf("SB2: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][y][(x+1)%Nx][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      printf("SB3: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      break;
      
   case 1:
      // upward staples
      stptmp  = L[t][z][(y+1)%Nx][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      printf("SF1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple = stptmp;
      stptmp  = L[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      printf("SF2: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = L[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      printf("SF3: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][z][(y+1)%Nx][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      printf("SB1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nx)%Nx][(y+1)%Nx][x][2]);
      stptmp *= conj(L[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= L[t][(z-1+Nx)%Nx][y][x][2];
      printf("SB2: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][(y+1)%Nx][x][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      printf("SB3: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      break;
      
   case 2:
      // upward staples
      stptmp  = L[t][(z+1)%Nx][y][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      printf("SF1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple = stptmp;
      stptmp  = L[t][(z+1)%Nx][y][x][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple += stptmp;
      printf("SF2: %f %f\n", creal(stptmp), cimag(stptmp));
      stptmp  = L[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      printf("SF3: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][(z+1)%Nx][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      printf("SB1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[t][(z+1)%Nx][(y-1+Nx)%Nx][x][1]);
      stptmp *= conj(L[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= L[t][z][(y-1+Nx)%Nx][x][1];
      printf("SB2: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][(z+1)%Nx][y][x][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      printf("SB3: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      break;
      
   case 3:
      // upward staples
      stptmp  = L[(t+1)%Nt][z][y][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      printf("SF1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple = stptmp;
      stptmp  = L[(t+1)%Nt][z][y][x][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      printf("SF2: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = L[(t+1)%Nt][z][y][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      printf("SF3: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[(t+1)%Nt][z][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      printf("SB1: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][z][(y-1+Nx)%Nx][x][1]);
      stptmp *= conj(L[t][z][(y-1+Nx)%Nx][x][dir]);
      stptmp *= L[t][z][(y-1+Nx)%Nx][x][1];
      printf("SB2: %f %f\n", creal(stptmp), cimag(stptmp));
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][(z-1+Nx)%Nx][y][x][2]);
      stptmp *= conj(L[t][(z-1+Nx)%Nx][y][x][dir]);
      stptmp *= L[t][(z-1+Nx)%Nx][y][x][2];
      printf("SB3: %f %f\n", creal(stptmp), cimag(stptmp));
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
   count = 0;
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {

      // p01, p02, p03
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][0]);
      stptmp *= conj(L[t][z][y][x][1]);
      p += stptmp;
      count++;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][0]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      count++;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][0]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      count++;

      // p12, p13
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][1]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      count++;
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][1]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      count++;

      // p23
      stptmp  = L[t][z][y][x][2];
      stptmp *= L[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][2]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      count++;
   }

	 }}}
   P = p;
   p /= Nx*Nx*Nx*Nt*6;
   printf("P %f %f, p %f %f\n", creal(P), cimag(P), creal(p), cimag(p));  
   printf("count %d\n", count);
   return creal(p);
}



int read_ascii_lat(char *fname) {
   char line[101];
   FILE *fp;
   int x,y,z,t,i,a,dir;
   float Ur, Ui;

   i = 0;
   fp = fopen(fname, "r");

   printf("reading from %s\n", fname);
   //   header
   while (i<3) {
      if( fgets (line, 100 , fp)!=NULL ) {
	 //	 fscanf(fp, "%s", &line);
	 i++;
	 //      printf("%s", i, line);
	 //puts(line);
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

}
		  
