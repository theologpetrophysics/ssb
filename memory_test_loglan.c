

#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <stdlib.h>



int allocateMemory(double** ptr, int n)
{
    *ptr = (double*)calloc(n, sizeof(double));  // allocate numFrame doubles
    if (*ptr == NULL) {
        fprintf(stderr, "calloc of size %d failed!\n", n);
        return 1;
    }
    else {
        fprintf(stderr, "\n Memory allocated  \n");
        return 0;
    }
}

