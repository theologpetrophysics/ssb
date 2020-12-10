

#include <stdio.h>
#include <string.h>
#include <math.h>  
#include <stdlib.h>



/* +++++++++++++++++++++++++++++++++++++++++++++
     Funtion to calculate average log value
+++++++++++++++++++++++++++++++++++++++++++++++*/
/*double smoothLogData(
    int numFrames,
    char* wtShape,
    double* logData)

{
    int i, lnumFrames;
    double sValue;
    double *wtFactor;
    double logDataSum;

    int debug = 1;

    lnumFrames = numFrames;

    if (debug == 1) {
        fprintf(stderr, "Number of frames is %d \n", numFrames);
        fprintf(stderr, "Weighting shape = %s \n", wtShape);
        fprintf(stderr, "Log data : ");
        for (i = 0; i < numFrames; i++) {
            fprintf(stderr, "%f, ", logData[i]);
        }
    }

    wtFactor = (double*)calloc(lnumFrames, sizeof(double));  // allocate numFrame doubles
    if (wtFactor == NULL) {
        fprintf(stderr, "calloc of size %d failed!\n", lnumFrames);   // could also call perror here
    }
    else {
        fprintf(stderr, "\n Memory allocated wtFactor  \n");
    }

    logDataSum = 0;

    for (i = 0; i < numFrames; i++) {
        fprintf(stderr,"i = %d \n", i);
       if (strncmp(wtShape, "box", 3) == 0) {
           fprintf(stderr,"1.\n");
           wtFactor[i] = 1.0;
           fprintf(stderr,"1.\n");
            fprintf(stderr, "wt factor %f, ", wtFactor[i]);
        } 


       logDataSum += logData[i] * wtFactor[i];
    }

    free(wtFactor);

    sValue = logDataSum / numFrames;
    if (debug == 1) {
        fprintf(stderr, "\n sValue = %f \n", sValue);
    }

    return sValue;

}*/

double smoothLogDataSimple(
    int numFrames,
    double* logData)

{
    int i;
    double sValue;
    double logDataSum;

    int debug = 0;

    if (debug == 1) {
        fprintf(stderr, "Number of frames is %d \n", numFrames);
        fprintf(stderr, "Log data : ");
        for (i = 0; i < numFrames; i++) {
            fprintf(stderr, "%f, ", logData[i]);
        }
    }

    logDataSum = 0;

    for (i = 0; i < numFrames; i++) {
        if (debug == 1) {
            fprintf(stderr, "i = %d \n", i);
        }

        logDataSum += logData[i];
    }


    sValue = logDataSum / numFrames;
    if (debug == 1) {
        fprintf(stderr, "\n sValue = %f \n", sValue);
    }

    return sValue;

}


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

