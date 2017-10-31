#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

complex get_staples(int x, int y, int z, int t, int dir);
void update();
double plaq();
double plaq2();
void printfC(char *s, complex z);
float magz(complex z);
void load_staples();



#define Nx 8
#define Nt 8

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

   if(argc != 5) {
      printf("Usage: %s beta warms trajecs measure\n", argv[0]);
      exit(0);
   } else {
      beta = atof(argv[1]);
      warms = atoi(argv[2]);
      trajecs = atoi(argv[3]);
      measure = atoi(argv[4]);
   }

   printf("# beta    = %f\n", beta);
   printf("# warms   = %d\n", warms);
   printf("# trajecs = %d\n", trajecs);
   printf("# measure = %d\n", measure);


   srand48(time(0));


   // cold lattice
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nx; z++) {
	 for(y=0; y<Nx; y++) {
	    for(x=0; x<Nx; x++) {
	       for(dir=0; dir<4; dir++) {
		  L[t][z][y][x][dir] = 1.0 + 0.0*I;
   }}}}}

   //#define BLAH
#ifdef BLAH
   c = cos(M_PI/3);
   s = sin(M_PI/3);

   L[0][0][7][1][1] = c - I*s;  // U(x+1)_y
   L[0][0][7][0][0] = c - I*s;  // U(x+1,y+1)_x
   L[0][0][0][0][1] = c + I*s;  // U(x)_y

   //   L[0][0][0][1][1] = 0 + I*0;  // U(x+1)_y
   //   L[0][0][1][0][0] = 0 - I*0;  // U(x+1,y+1)_x
   //   L[0][0][0][0][1] = 0 - I*0;  // U(x)_y

   stp = L[0][0][0][1][1] * conj(L[0][0][1][0][0]) * conj(L[0][0][0][0][1]);
   printfC("explicit stp=", stp);

   stp = get_staples(0,0,0,0,0);
   printfC("get_stp=", stp);
   exit(0);
#endif


   //   plaquette = plaq2();
   //printf("PLAQ %f %f\n", beta, plaquette);

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
   complex staple, stptmp, p;
   int x,y,z,t,dir,dir1;

   p = 0.0 + 0.0*I;
   
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
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][0]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][0]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      
      // p12, p13
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][2];
      stptmp *= conj(L[t][(z+1)%Nx][y][x][1]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Nx][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][1]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;

      // p23
      stptmp  = L[t][z][y][x][2];
      stptmp *= L[t][(z+1)%Nx][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][2]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
   }

	 }}}

   p /= Nx*Nx*Nx*Nt*6;

   return creal(p);
}


