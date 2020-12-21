

#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <stdlib.h>

#include "theolog_general_functions.h"

#define LGROUPSIZE 5000

/*
struct inputLogData {

    double* depth;
    double* tvd;
    double* gr;
    double* vsh;
};

*/



void ssbFromLoglan(
    int numpts,
    double depthAveWindow,
    double *depthLog,
    double *tvdLog,
    double *grLog, 
    double *vshLog,
    double lithGroupMinThick,
    double *vshSmth,
    double* vshFirstDeriv,
    double* vshSecondDeriv,
    double* vshSmthFirstDeriv,
    double* vshSmthSecondDeriv,
    double * lithLogValue)

{
    int i, k, j, grpcnt;

    int halfSmthWindow;
    double sampleRate;
    double* lithGroupDepth;
    double* lithGroupThick;
    double* lithGroupVsh;
    double* lithGroupValue;
    int numLithGroups = 0;


    allocateMemory1DD(&lithGroupDepth, numpts, 0);
    allocateMemory1DD(&lithGroupThick, numpts,0);
    allocateMemory1DD(&lithGroupVsh, numpts, 0);
    allocateMemory1DD(&lithGroupValue, numpts, 0);

    char lithGroupName[LGROUPSIZE][10] = { {0} };
    char lithGroupLith[LGROUPSIZE][5] = { {0} };

   /* struct inputLogData inData;
    inData.depth = depthLog;
    inData.tvd = tvdLog;
    inData.gr = grLog;
    inData.vsh = vshLog;
    */

    fprintf(stderr, "processing range (%f - %f m) \n", depthLog[0], depthLog[numpts-1]);


    // calculate number of frames based on input data 
    sampleRate = (depthLog[numpts - 1] - depthLog[0]) / numpts;

    // set half smoothing window
    halfSmthWindow = (int)fmin( (depthAveWindow / 2.0 ) / sampleRate, 499);

    // inform the user ******************
    fprintf(stderr, "Sample rate of data is %f.3 (m) \n", sampleRate);
    fprintf(stderr, "Number of log samples to process is : %i \n The smoothfactor halfwindow is : %i \n", numpts, halfSmthWindow);
    //***********************************

    // check that Vsh is limited (0<=1)
    for (i = 0; i < numpts; i++) {
        vshLog[i] = fmin(vshLog[i], 1.0);
        vshLog[i] = fmax(vshLog[i], 0.0);
    }



    //call main SSB workflow
    ssbMain(
        numpts,
        halfSmthWindow,
        depthLog,
        tvdLog,
        grLog,
        vshLog,
        vshSmth,
        vshFirstDeriv,
        vshSecondDeriv,
        vshSmthFirstDeriv,
        vshSmthSecondDeriv,
        lithLogValue,
        lithGroupMinThick,
        lithGroupValue,
        lithGroupDepth,
        lithGroupThick,
        lithGroupVsh,
        &numLithGroups);

    for (i = 0; i < numLithGroups-1; i++) {
        if (lithGroupValue[i] == 1) {
            //fprintf(stderr, "lithgroup %d - ", i);
            strncpy(lithGroupLith[i], "sh", 2);
        }
        else {
            strncpy(lithGroupLith[i], "ss", 2);
        }       
    }


    grpcnt = 1;
    for (i = numLithGroups - 2; i >= 0; i--, grpcnt++) {
        sprintf(lithGroupName[i], "LG-%d", grpcnt);
    }

    //write out results to .csv file
    FILE* vshlgfile;

    //fprintf(stderr, "numgrps = %d \n", *numLithGroups);

    vshlgfile = fopen("./data/vshlithogroup.csv", "w");
    if (vshlgfile != NULL) {
        fprintf(stderr, "./data/vshlithogroup.csv has been opened for write \n");
    }

    fprintf(stderr, "numgrps = %d \n", numLithGroups);

    for (i = 0; i < numLithGroups; i++) {
        fprintf(vshlgfile, "%f, %f, %s, %s, %f \n", lithGroupDepth[i], lithGroupThick[i], lithGroupLith[i], lithGroupName[i], lithGroupVsh[i]);
    }

    fclose(vshlgfile);


    // return to loglan


}
