

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
    char* lithgroupMethod,
    char* ecsMethod,
    double ecsShaleBreakLimit,
    char* tsfPickMethod,
    char* optElementLog,
    double lithGroupMinThick,
    double* vshSmth,
    double* grSmth,
    double* firstDeriv,
    double* secondDeriv,
    double* lithLogValue)

{
    int i, k, j, grpcnt;

    int halfSmthWindow;
    double sampleRate;

    //int numLithGroups = 0;


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
    //halfSmthWindow = fmax((double)halfSmthWindow, 1.0);

    // inform the user ******************
    fprintf(stderr, "Sample rate of data is %f (m) \n", sampleRate);
    //***********************************

    // check that Vsh is limited (0<=1)
    for (i = 0; i < numpts; i++) {
        vshLog[i] = fmin(vshLog[i], 1.0);
        vshLog[i] = fmax(vshLog[i], 0.0);
    }

    fprintf(stderr, "before deriv... \n \n");

    //call main SSB workflow
    ssbMain(
        numpts,
        halfSmthWindow,
        depthLog,
        tvdLog,
        grLog,
        vshLog,
        lithgroupMethod,
        ecsMethod,
        ecsShaleBreakLimit,
        tsfPickMethod,
        optElementLog,
        vshSmth,
        grSmth,
        firstDeriv,
        secondDeriv,
        lithLogValue,
        lithGroupMinThick);



}
