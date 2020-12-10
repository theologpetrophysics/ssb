

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "theolog_general_functions.h"



void ssbMain(
    int numpts,
    int halfSmthWindow,
    double* depth,
    double* tvd,
    double* gr,
    double* vsh,
    double* vshSmth )

{
    int i, k, j;
    double smthData[1000] = { 0 }, localsValue;
    int framesToBeSmth;

    int debug = 0;

    if (debug == 1) {
        // inform the user ******************
        fprintf(stderr, "2. Number of log samples to process is : %i \n 2. The smoothfactor halfwindow is : %i \n", numpts, halfSmthWindow);
        //***********************************

        // Calculate smoothed VSH curve *****
        fprintf(stderr, "generating smoothed Vsh curve");
        //***********************************
    }

    for (i = 0; i < numpts; i++) {
        j = 0;
        //generate array of data to be averaged
        for (k = fmax(i- halfSmthWindow,0); k <= fmin((i + halfSmthWindow),numpts-1); k++, j++) {
            smthData[j] = vsh[k];
            if (debug == 1) {
                fprintf(stderr, "%f, ", smthData[j]);
            }
         }
        framesToBeSmth = j;
        if (debug == 1) {
            fprintf(stderr, "frames to be smoothed %d; \n", framesToBeSmth);
        }


        // call smoothing function

        localsValue = smoothLogDataSimple(framesToBeSmth, smthData);
        vshSmth[i] = localsValue;
        if (debug == 1) {
            fprintf(stderr, "Smoothed value = %f \n ", localsValue);
        }
    }


}
