

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "theolog_general_functions.h"

#define MISSING -999.25

void ssbMain(
    int numpts,
    int halfSmthWindow,
    double* depth,
    double* tvd,
    double* gr,
    double* vsh,
    char* lithgoupMethod,
    double* vshSmth,
    double* vshFirstDeriv,
    double* vshSecondDeriv,
    double* vshSmthFirstDeriv,
    double* vshSmthSecondDeriv,
    double* lithValue,
    double lithGroupMinThick,
    double* lithGroupValue,
    double* lithGroupDepth,
    double* lithGroupThick,
    double* lithGroupVsh,
    int* numLithGroups)

{
    int i, k, j, lithframecnt;
    double smthData[1000] = { 0 }, localsValue;
    int framesToBeSmth;

    double liththicksum, liththickmin, liththickmax, lithvshsum;

    int debug = 0;

    //allocateMemory1DI(&lithValue, numpts, 0);

    // inform the user ******************
    fprintf(stderr, "2. Number of log samples to process is : %i \n 2. The smoothfactor halfwindow is : %i \n", numpts, halfSmthWindow);
    //***********************************

    // Calculate smoothed VSH curve *****
    fprintf(stderr, "generating smoothed Vsh curve \n");
    //***********************************

    // logic to smooth the VSH curve
    //Cycle through data frames
    for (i = 0; i < numpts; i++) {
        j = 0;
        //generate array of data to be averaged
        for (k = fmax(i - halfSmthWindow, 0); k <= fmin((i + halfSmthWindow), numpts - 1); k++, j++) {
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

        localsValue = smoothLogData(framesToBeSmth, "box", smthData);
        vshSmth[i] = localsValue;
        if (debug == 1) {
            fprintf(stderr, "Smoothed value = %f \n ", localsValue);
        }
    }

    // caluclate first and second derivative of Vsh and smoothed Vsh
    for (i = 0; i < numpts - 1; i++) {
        vshFirstDeriv[i] = (vsh[i + 1] - vsh[i]) / (depth[i + 1] - depth[i]);
        vshSmthFirstDeriv[i] = (vshSmth[i + 1] - vshSmth[i]) / (depth[i + 1] - depth[i]);
    }

    for (i = 0; i < numpts - 1; i++) {
        vshSecondDeriv[i] = (vshFirstDeriv[i + 1] - vshFirstDeriv[i]) / (depth[i + 1] - depth[i]);
        vshSmthSecondDeriv[i] = (vshSmthFirstDeriv[i + 1] - vshSmthFirstDeriv[i]) / (depth[i + 1] - depth[i]);
    }


    // logic to calculate vsh lithology groups
    // First pick sand/shale depths
    for (i = 0; i < numpts; i++) {

        if (vsh[i] == vshSmth[i] == 1.0) {
            lithValue[i] = 1;
        }
        else if (strncmp(lithgoupMethod, "AVERAGE", 7) == 0) {
            if (vsh[i] >= vshSmth[i]) {
                //shale lithogroup
                lithValue[i] = 1;
            }
            else {
                // sand lithogroup
                lithValue[i] = 0;
            }
        }
        else {
            if (vshSecondDeriv[i] >= 0) {
                lithValue[i] = 0;
            }
            else {
                lithValue[i] = 1;
            }
        }
    }

    //lithgroup counter
    j = 0;
    liththicksum = 0;
    liththickmin = 100;
    liththickmax = 0;
    lithvshsum = 0;
    lithframecnt = 0;

    // start with top frame and define group tops
    lithGroupDepth[j] = depth[0];
    if (lithValue[0] == 1) {
        lithGroupValue[j] = 1;
    }
    else {
        lithGroupValue[j] = 0;
    }

    lithvshsum += vsh[0];
    lithframecnt++;

    for (i = 1; i < numpts-1; i++) {
        if (lithValue[i] != lithValue[i - 1] &
            lithValue[i] == lithValue[i+1] &
            lithValue[i-1] == lithValue[i-2] &
            (depth[i] - lithGroupDepth[j]) >= lithGroupMinThick ) {

            // lithgroup bed statistics of previous bed
            lithGroupThick[j] = depth[i] - lithGroupDepth[j];
            liththicksum += lithGroupThick[j];
            if (lithGroupThick[j] > liththickmax) liththickmax = lithGroupThick[j];
            if (lithGroupThick[j] < liththickmin) liththickmin = lithGroupThick[j];

            lithGroupVsh[j] = lithvshsum / lithframecnt;

            // no increment and define current bed lith
            j++;
            lithGroupValue[j] = lithValue[i];
            lithGroupDepth[j] = depth[i];

            if (debug == 1) {
                fprintf(stderr, "grp depth  = %f \n ", lithGroupDepth[j]);
            }

            lithframecnt = 1;
            lithvshsum = 0;

        }
        else {
            lithvshsum += vsh[i];
            lithframecnt++;
        }
    }

    i = numpts-1;
    // last frame of log data
    lithGroupThick[j] = depth[i] - lithGroupDepth[j];
    liththicksum += lithGroupThick[j];
    if (lithGroupThick[j] > liththickmax) liththickmax = lithGroupThick[j];
    if (lithGroupThick[j] < liththickmin) liththickmin = lithGroupThick[j];
    lithGroupValue[j] = lithValue[i - 1];
    lithGroupVsh[j] = lithvshsum / lithframecnt;
    j++;
    lithGroupDepth[j] = depth[i];
    lithGroupThick[j] = MISSING;
    lithGroupVsh[j] = MISSING;
    fprintf(stderr, "\n **** j = %d, grp depth  = %f  ***\n ", j, lithGroupDepth[j]);
    lithframecnt = 0;
    lithvshsum = 0;

   
    //j++;
    //lithGroupDepth[j] = depth[numpts-1];
    if (debug == 1) {
        fprintf(stderr, "grp base depth  = %f \n ", lithGroupDepth[j]);
    }


    
    *numLithGroups = j + 1;
    fprintf(stderr, "Number of Vsh Lithology Groups = %d \n", *numLithGroups-1);
    fprintf(stderr, "Average thickness of Vsh Lithology Groups = %f \n", liththicksum / *numLithGroups);
    fprintf(stderr, "Minimum thickness of Vsh Lithology Groups = %f \n", liththickmin);
    fprintf(stderr, "Maximum thickness of Vsh Lithology Groups = %f \n", liththickmax);


}
