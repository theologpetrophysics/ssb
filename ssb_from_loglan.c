

#include <stdio.h>
#include <string.h>
#include <math.h>   


void ssbFromLoglan(
    int numpts,
    int smthFactor,
    double *depthLog,
    double *tvdLog,
    double *grLog, 
    double *vshLog,
    double *vshSmth )

{
    int i, k, j;
    double smthData[1000];
    int halfSmthWindow, framesToBeSmth;

    // set this initally - this is half smoothing window
    halfSmthWindow = smthFactor;

    /*
    // find number of depths in half smthFactor
    if (smthFactor % 2 == 0) {
        framesToBeSmth = smthFactor / 2;
    }
    else {
        framesToBeSmth = (smthFactor -1) / 2;
    }
    */


    fprintf(stderr, "Number of points is : %i, The smoothfactor is : %i \n", numpts, smthFactor);

    // test to display for 5 data frames
    for (i = 0; i < 5; i++) {
        fprintf(stderr, "depth = %f \n", depthLog[i]);
        fprintf(stderr, "tvd = %f \n", tvdLog[i]);
        fprintf(stderr, "gr = %f \n", grLog[i]);
        fprintf(stderr, "vsh = %f \n", vshLog[i]);

        j = 0;
        //generate array of data to be averaged
        for (k = fmax(i- halfSmthWindow,0); k <= (i + halfSmthWindow); k++, j++) {
            smthData[j] = vshLog[k];
            fprintf(stderr, "data for smoothing = %f \n", smthData[j]);
         }
        framesToBeSmth = j;
        fprintf(stderr, "frames to be smoothed = %i \n", framesToBeSmth);

    }





}


double smoothLogData(
    int numFrames,
    char* wtShape,
    double* logData,
    double sValue)

{
    int i;
    //double *wtFactor[numFrames];

    for (i = 0; i < numFrames; i++) {
       /* if (strncmp(wtShape, "box", 3) == 0) {
            wtFactor[i] = 1;
        } */

        fprintf(stderr, "%f \n", logData[i]);
    }

    return sValue;

}

/*
int fmin(int x, int y)
{
    return (x < y) ? x : y;
}

int fmax(int x, int y)
{
    return (x > y) ? y : x;
}

*/
