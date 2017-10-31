////////////////////////////////////////////////
// U(1) Measurements
////////////////////////////////////////////////

#include "u1includes.h"



double plaq_twist() {
   complex staple, stptmp, P, p;
   int x,y,z,t,dir,dir1;
   int count;
   complex P01,P02,P03,P12,P13,P23;
   complex wphase, wphasecc;


   p = 0.0 + 0.0*I;
   // count = 0;
   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {

      // put in BC phases for U_x = L[..][0]
      wphase = (y==Ny-1) ? cos(Ny*w) - I*sin(Ny*w) : 1.0 + I*0.0;
      wphasecc = (y==0) ? cos(Ny*w) + I*sin(Ny*w) : 1.0 + I*0.0;

      // p01, p02, p03
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj(wphase * L[t][z][(y+1)%Ny][x][0]);  // w
      stptmp *= conj(L[t][z][y][x][1]);
      p += stptmp;
      P01 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][0];
      stptmp *= L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][0]);
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
      stptmp *= L[t][z][(y+1)%Ny][x][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][1]);
      stptmp *= conj(L[t][z][y][x][2]);
      p += stptmp;
      P12 = stptmp;
      //    count++;
      stptmp  = L[t][z][y][x][1];
      stptmp *= L[t][z][(y+1)%Ny][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][1]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P13 = stptmp;
      //    count++;

      // p23
      stptmp  = L[t][z][y][x][2];
      stptmp *= L[t][(z+1)%Nz][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][2]);
      stptmp *= conj(L[t][z][y][x][3]);
      p += stptmp;
      P23 = stptmp;
      //    count++;

      /**
      printf("P %d %d %d %d :\n(%f,%f)\n(%f,%f)\n(%f,%f)\n(%f,%f)\n(%f,%f)\n(%f,%f)\n",
	     x,y,z,t, 
	     creal(P01),cimag(P01), creal(P02),cimag(P02),
	     creal(P12),cimag(P12), creal(P03),cimag(P03), 
	     creal(P13),cimag(P13), creal(P23),cimag(P23));
      **/
   }}}}

   p /= Nvol*6;
   //   printf("P %f %f, p %f %f\n", creal(P), cimag(P), creal(p), cimag(p));  
   //   printf("count %d\n", count);
   return creal(p);
}




double plaq() {
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



double plaqfromstaple() {
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


