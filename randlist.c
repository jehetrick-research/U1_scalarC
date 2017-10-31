#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main() {
   int n=1000000;
   int i;

   srand48(1000);

   
   for(i=0; i<n; i++) {
      printf("%f\n", drand48());
   }

   return 0;
}
