//////////////////////////////////////////////
// Compute Staples
//////////////////////////////////////////////

#include "u1includes.h"



void load_staples_twist(int dir) {
   int x,y,z,t;

   for(t=0; t<Nt; t++) {
      for(z=0; z<Nz; z++) {
	 for(y=0; y<Ny; y++) {
	    for(x=0; x<Nx; x++) {
		  S[t][z][y][x][dir] = get_staples_twist(x,y,z,t,dir);
   }}}}
}


complex get_staples_twist(int x, int y, int z, int t, int dir) {
   complex staple, stptmp;
   complex wphase, wphasecc;

   // put in BC phases for U_x = L[..][0]
   wphase = (y==Ny-1) ? cos(Ny*w) - I*sin(Ny*w) : 1.0 + I*0.0;
   wphasecc = (y==0) ? cos(Ny*w) + I*sin(Ny*w) : 1.0 + I*0.0;

   switch(dir) {
   case 0:
      // upward staples
      stptmp  = L[t][z][y][(x+1)%Nx][1];
      stptmp *= conj( wphase * L[t][z][(y+1)%Ny][x][0]);  // w
      stptmp *= conj(L[t][z][y][x][1]);
      staple = stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][0]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][0]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][z][(y-1+Ny)%Ny][(x+1)%Nx][1]);
      stptmp *= conj(wphasecc * L[t][z][(y-1+Ny)%Ny][x][0]); // w*
      stptmp *= L[t][z][(y-1+Ny)%Ny][x][1];
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nz)%Nz][y][(x+1)%Nx][2]);
      stptmp *= conj(L[t][(z-1+Nz)%Nz][y][x][0]);
      stptmp *= L[t][(z-1+Nz)%Nz][y][x][2];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][y][(x+1)%Nx][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][0]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 1:
      // upward staples
      stptmp  = wphase * L[t][z][(y+1)%Ny][x][0];  // w
      stptmp *= conj(L[t][z][y][(x+1)%Nx][1]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[t][z][(y+1)%Ny][x][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][1]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = L[t][z][(y+1)%Ny][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][1]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(wphase * L[t][z][(y+1)%Ny][(x-1+Nx)%Nx][0]); // w
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][1]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nz)%Nz][(y+1)%Ny][x][2]);
      stptmp *= conj(L[t][(z-1+Nz)%Nz][y][x][1]);
      stptmp *= L[t][(z-1+Nz)%Nz][y][x][2];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][(y+1)%Ny][x][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][1]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 2:
      // upward staples
      stptmp  = L[t][(z+1)%Nz][y][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][2]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[t][(z+1)%Nz][y][x][1];
      stptmp *= conj(L[t][z][(y+1)%Nx][x][2]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = L[t][(z+1)%Nz][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][2]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][(z+1)%Nz][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][2]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[t][(z+1)%Nz][(y-1+Ny)%Ny][x][1]);
      stptmp *= conj(L[t][z][(y-1+Ny)%Ny][x][2]);
      stptmp *= L[t][z][(y-1+Ny)%Ny][x][1];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][(z+1)%Nz][y][x][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][2]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 3:
      // upward staples
      stptmp  = L[(t+1)%Nt][z][y][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][3]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[(t+1)%Nt][z][y][x][1];
      stptmp *= conj(L[t][z][(y+1)%Ny][x][3]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = L[(t+1)%Nt][z][y][x][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][3]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[(t+1)%Nt][z][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][3]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][z][(y-1+Ny)%Ny][x][1]);
      stptmp *= conj(L[t][z][(y-1+Ny)%Ny][x][3]);
      stptmp *= L[t][z][(y-1+Ny)%Ny][x][1];
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][(z-1+Nz)%Nz][y][x][2]);
      stptmp *= conj(L[t][(z-1+Nz)%Nz][y][x][3]);
      stptmp *= L[t][(z-1+Nz)%Nz][y][x][2];
      staple += stptmp;
      break;
   }
    
   return staple;
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
      stptmp *= conj(L[t][z][(y+1)%Ny][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple = stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = L[t][z][y][(x+1)%Nx][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][z][(y-1+Ny)%Ny][(x+1)%Nx][1]);
      stptmp *= conj(L[t][z][(y-1+Ny)%Ny][x][dir]);
      stptmp *= L[t][z][(y-1+Ny)%Ny][x][1];
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nz)%Nz][y][(x+1)%Nx][2]);
      stptmp *= conj(L[t][(z-1+Nz)%Nz][y][x][dir]);
      stptmp *= L[t][(z-1+Nz)%Nz][y][x][2];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][y][(x+1)%Nx][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 1:
      // upward staples
      stptmp  = L[t][z][(y+1)%Ny][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[t][z][(y+1)%Ny][x][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      stptmp  = L[t][z][(y+1)%Ny][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][z][(y+1)%Ny][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[t][(z-1+Nz)%Nz][(y+1)%Ny][x][2]);
      stptmp *= conj(L[t][(z-1+Nz)%Nz][y][x][dir]);
      stptmp *= L[t][(z-1+Nz)%Nz][y][x][2];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][z][(y+1)%Ny][x][3]);
      stptmp *= conj(L[(t-1+Nt)%Nt][z][y][x][dir]);
      stptmp *= L[(t-1+Nt)%Nt][z][y][x][3];
      staple += stptmp;
      break;
      
   case 2:
      // upward staples
      stptmp  = L[t][(z+1)%Nz][y][x][0];
      stptmp *= conj(L[t][z][y][(x+1)%Nx][dir]);
      stptmp *= conj(L[t][z][y][x][0]);
      staple = stptmp;
      stptmp  = L[t][(z+1)%Nz][y][x][1];
      stptmp *= conj(L[t][z][(y+1)%Ny][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = L[t][(z+1)%Nz][y][x][3];
      stptmp *= conj(L[(t+1)%Nt][z][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][3]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[t][(z+1)%Nz][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[t][(z+1)%Nz][(y-1+Ny)%Ny][x][1]);
      stptmp *= conj(L[t][z][(y-1+Ny)%Ny][x][dir]);
      stptmp *= L[t][z][(y-1+Ny)%Ny][x][1];
      staple += stptmp;
      stptmp  = conj(L[(t-1+Nt)%Nt][(z+1)%Nz][y][x][3]);
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
      stptmp *= conj(L[t][z][(y+1)%Ny][x][dir]);
      stptmp *= conj(L[t][z][y][x][1]);
      staple += stptmp;
      stptmp  = L[(t+1)%Nt][z][y][x][2];
      stptmp *= conj(L[t][(z+1)%Nz][y][x][dir]);
      stptmp *= conj(L[t][z][y][x][2]);
      staple += stptmp;
      
      // downward staples
      stptmp  = conj(L[(t+1)%Nt][z][y][(x-1+Nx)%Nx][0]);
      stptmp *= conj(L[t][z][y][(x-1+Nx)%Nx][dir]);
      stptmp *= L[t][z][y][(x-1+Nx)%Nx][0];
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][z][(y-1+Ny)%Ny][x][1]);
      stptmp *= conj(L[t][z][(y-1+Ny)%Ny][x][dir]);
      stptmp *= L[t][z][(y-1+Ny)%Ny][x][1];
      staple += stptmp;
      stptmp  = conj(L[(t+1)%Nt][(z-1+Nz)%Nz][y][x][2]);
      stptmp *= conj(L[t][(z-1+Nz)%Nz][y][x][dir]);
      stptmp *= L[t][(z-1+Nz)%Nz][y][x][2];
      staple += stptmp;
      break;
   }
    
   return staple;
}
