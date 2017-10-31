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
