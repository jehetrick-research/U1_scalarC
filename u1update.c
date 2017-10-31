/////////////////////////////////////////
// U(1) Update code
/////////////////////////////////////////

#include "u1includes.h"



///////////////////////////////////////////////////////////
// Do NOT update dir[0] links on xy plane @z=z_xyplane 
///////////////////////////////////////////////////////////

void update_flux(int z_xyplane) {
   int x,y,z,t,dir, i;
   double r, th, scale, oldaction, newaction, actiondiff;
   complex staples, U, Unew, dU;
   double latmag, accrej, pold, pnew, pdiff, accinc, accyes;

   scale = 1.4;

   latmag = 0.0;
   accrej = 0.0;


   for(dir=0; dir<4; dir++) {
      load_staples_twist(dir);

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
      load_staples_twist(dir);

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




