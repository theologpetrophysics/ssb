/* This is an example of a simple call to external c

    call it as cvshgr ( gr, gr_ma, gr_sh, &vsh_gr )

*/

#include <stdio.h>


void cvshgr ( double *gr, 
              double gr_ma, 
              double gr_sh, 
              double *vsh_gr 
            )
{


    *vsh_gr = ( gr - gr_ma ) / ( gr_sh - gr_ma );
    fprintf(stderr, "vsh = %f \n", *vsh_gr);
  

  return;
}
