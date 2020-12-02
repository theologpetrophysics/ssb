/* This is an example of a simple call to external c

    call it as cvshgr ( gr, gr_ma, gr_sh, &vsh_gr )

*/

#include <stdio.h>


void cvshgr ( int numpts,
              double *gr, 
              double gr_ma, 
              double gr_sh,
              double *vshGr )
{


int i;

for (i=0;i<numpts; i++) {
   vshGr[i] = ( gr[i] - gr_ma ) / ( gr_sh - gr_ma );
    //fprintf(stderr, "vsh = %f \n", vshGr[i]);
}

   
  

  return;
}
