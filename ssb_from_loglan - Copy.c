

#include <stdio.h>
#include <string.h>
#include <math.h>  


double smoothLogData(
    int numFrames,
    char* wtShape,
    double* logData);




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
    double smthData[1000] = { 0 };
    int halfSmthWindow, framesToBeSmth;
    double sampleRate;

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

    // Calculate smoothed VSH curve *****
    fprintf(stderr, "generating smoothed Vsh curve");
    //***********************************
    for (i = 0; i < numpts; i++) {
        j = 0;
        //generate array of data to be averaged
        for (k = fmax(i- halfSmthWindow,0); k <= fmin((i + halfSmthWindow),numpts-1); k++, j++) {
            smthData[j] = vshLog[k];

         }
        framesToBeSmth = j;

        // call smoothing function
        vshSmth[i] = smoothLogData(framesToBeSmth, "box", smthData);
    }


}


double smoothLogData(
    int numFrames,
    char* wtShape,
    double* logData)

{
    int i;
    double sValue;
    double wtFactor[1000] = { 0 };
    double logDataSum;

    logDataSum = 0;

    for (i = 0; i < numFrames; i++) {
       if (strncmp(wtShape, "box", 3) == 0) {
            wtFactor[i] = 1.0;
        } 

       logDataSum += logData[i] * wtFactor[i];
    }

    sValue = logDataSum / numFrames;

    return sValue;

}




/* unused code

        //fprintf(stderr, "depth = %f \n", depthLog[i]);
        //fprintf(stderr, "tvd = %f \n", tvdLog[i]);
        //fprintf(stderr, "gr = %f \n", grLog[i]);
        //fprintf(stderr, "vsh = %f \n", vshLog[i]);

        //fprintf(stderr, "frames to be smoothed = %i \n", framesToBeSmth);
        //fprintf(stderr, "data for smoothing = %f \n", smthData[j]);


*/

