

#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <stdlib.h>

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
    double *vshSmth )

{
    int i, k, j;

    int halfSmthWindow;
    double sampleRate;

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
        vshSmth);

    // return to loglan


}
