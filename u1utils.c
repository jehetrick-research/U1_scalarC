/////////////////////////////////////////////////////
// U(1) Functions
/////////////////////////////////////////////////////

#include "u1includes.h"




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

void dumpplaq_twist() {
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
	    // put in BC phases for U_x = L[..][0]
	    wphase = (y==Ny-1) ? cos(Ny*w) - I*sin(Ny*w) : 1.0 + I*0.0;
	    wphasecc = (y==0) ? cos(Ny*w) + I*sin(Ny*w) : 1.0 + I*0.0;

	    for(x=0; x<Nx; x++) {

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

      /*
      printf("%d %d %d %d : (%f,%f) (%f,%f) (%f,%f)\n\t\t(%f,%f) (%f,%f) (%f,%f)\n",
	     x,y,z,t, 
	     creal(P01),cimag(P01), creal(P02),cimag(P02),
	     creal(P03),cimag(P03), creal(P12),cimag(P12),
	     creal(P13),cimag(P13), creal(P23),cimag(P23));
      */

      printf(
"%d %d %d %d P01 %f %f %f P02 %f %f P03 %f %f P12 %f %f P13 %f %f P23 %f %f\n",
	     x,y,z,t, 
	     creal(P01),cimag(P01), atan2(cimag(P01),creal(P01)),
	     creal(P02),cimag(P02),
	     creal(P03),cimag(P03), creal(P12),cimag(P12),
	     creal(P13),cimag(P13), creal(P23),cimag(P23));
	     

   }}}}
}



////////////////////////////////////////////////////
// Read an ASCII Milc Lattice
////////////////////////////////////////////////////


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

   fclose(fp);

}
		  
