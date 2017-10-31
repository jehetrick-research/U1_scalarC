#ifndef _U1DEFS
#define _U1DEFS

#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef struct {
   int x;
   int y;
   int z;
   int t;
} parity;


#define EXTERN extern

// Global Variables
///////////////////
EXTERN complex *****L;
EXTERN complex *****S;
EXTERN double beta;
EXTERN int Nvol;
EXTERN int Np; // number of EVEN, ODD sites: 1/2 Nvol
EXTERN int Nx; 
EXTERN int Ny;
EXTERN int Nz;
EXTERN int Ns; // N_space
EXTERN int Nt; //N_time
EXTERN parity *EVEN;
EXTERN parity *ODD;
EXTERN double **mflux;
EXTERN double w; // twist angle




complex get_staples(int x, int y, int z, int t, int dir);
complex get_staples_twist(int x, int y, int z, int t, int dir);
void update();
double plaq();
double plaqfromstaple();
double plaq_twist();
void printfC(char *s, complex z);
float magz(complex z);
void load_staples(int dir);
void load_staples_twist(int dir);
int read_ascii_lat(char *fname);
void fluxplane(int z_xyplane);
void fluxplane_twisttwist(int z_xyplane);
void update_flux(int z_xyplane);
void update_flux_twist(int z_xyplane);
void meas_flux();
void meas_flux_twist();
void dumplat();
void dumpplaq();
void dumpplaq_twist();

#define COLD 0
#define HOT 1
#define LOADFILE 2





#endif
