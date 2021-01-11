

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "theolog_general_functions.h"

#define MISSING -999.25
#define LGROUPSIZE 5000


void charArrayPrintTest(
    int n,
    char a[][10]
)
{

    for (int i = 0; i < n; i++) {
        fprintf(stderr, "test: %s \n ", a[i]);
    }

}


void ssbMain(
    int numpts,
    int halfSmthWindow,
    double* depth,
    double* tvd,
    double* gr,
    double* vsh,
    char* lithGroupMethod,
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
    int i, k, j, lithframecnt, grpcnt;
    int numECS;
    double smthData[1000] = { 0 }, localsValue;
    int framesToBeSmth;

    double liththicksum, liththickmin, liththickmax, lithvshsum;

    double* ecsDepth;

    int debug = 0;

    //allocateMemory1DI(&lithValue, numpts, 0);

    char lithGroupName[LGROUPSIZE][10] = { {0} };
    char lithGroupLith[LGROUPSIZE][5] = { {0} };
    char ecColour[LGROUPSIZE][10] = { {0} };
    char ecsName[LGROUPSIZE][10] = { {0} };
    char ecsColour[LGROUPSIZE][10] = { {0} };

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

    //-----------------------------------------------------------------------------------------------------------------------------
    // logic to calculate element sets (put this in function)----------------------------------------------------------------------

    // First pick sand/shale depths
    for (i = 0; i < numpts; i++) {

        if (vsh[i] == vshSmth[i] == 1.0) {
            lithValue[i] = 1;
        }
        else if (strncmp(lithGroupMethod, "AVERAGE", 7) == 0) {
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

            lithframecnt = 0;
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

    if (debug == 1) {
        fprintf(stderr, "grp base depth  = %f \n ", lithGroupDepth[j]);
    }
   
    *numLithGroups = j + 1;
    fprintf(stderr, "Number of Element Sets = %d \n", *numLithGroups-1);
    fprintf(stderr, "Average thickness of Element Sets = %f \n", liththicksum / *numLithGroups);
    fprintf(stderr, "Minimum thickness of Element Sets = %f \n", liththickmin);
    fprintf(stderr, "Maximum thickness of Element Sets = %f \n", liththickmax);

    // generate lithology of Element sets
    for (i = 0; i < *numLithGroups - 1; i++) {
        if (lithGroupValue[i] == 1) {
            //fprintf(stderr, "lithgroup %d - ", i);
            strncpy(lithGroupLith[i], "sh", 2);
            strncpy(ecColour[i], "gray15", 6);
        }
        else {
            strncpy(lithGroupLith[i], "ss", 2);
            strncpy(ecColour[i], "yellow", 6);
        }
    }
       
    // generate name of element sets
    grpcnt = 1;
    for (i = *numLithGroups - 2; i >= 0; i--, grpcnt++) {
        sprintf(lithGroupName[i], "ES-%d", grpcnt);        
    }
  
    //write out results to .csv file
    FILE* vshlgfile;

    //fprintf(stderr, "numgrps = %d \n", *numLithGroups);

    vshlgfile = fopen("./data/vshlithogroup.csv", "w");
    if (vshlgfile != NULL) {
        fprintf(stderr, "./data/vshlithogroup.csv has been opened for write \n");
    }

    fprintf(stderr, "numgrps = %d \n", *numLithGroups);

    grpcnt = 1;
    for (i = *numLithGroups - 2; i >= 0; i--) {
        fprintf(stderr, "name = %s \n", lithGroupName[i]);
    }

    for (i = 0; i < *numLithGroups; i++) {
        fprintf(vshlgfile, "%f, %f, %s, %s, %s, %f \n", lithGroupDepth[i], lithGroupThick[i], lithGroupLith[i], lithGroupName[i], ecColour[i], lithGroupVsh[i]);
    }

    fclose(vshlgfile);

    //  end of potential sub function for generating the element set data----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------



    //  Now need to find the shale breaks....
    //  start from base, find when shale is thicker above than below
    //  is shale above current surface is thicker than next shale belo then current surface is a shale break
    //  Zone between shale breaks is called Element complex sets (ECS)

    // top surface is treated as top of ECS

    //lithgroup counter

    allocateMemory1DD(&ecsDepth, *numLithGroups, 0);


    k = 0;
   /* ecsthicksum = 0;
    ecsthickmin = 100;
    ecsthickmax = 0;
    ecsvshsum = 0;
    ecframecnt = 0;*/


   for (j = 0; j < *numLithGroups; j++) {
       if (j == 0) {
           ecsDepth[k] = lithGroupDepth[j];
           fprintf(stderr, "Top at  = %f \n ", ecsDepth[k]);
           k++;
       }
       else if (j == *numLithGroups - 1) {
           ecsDepth[k] = lithGroupDepth[j];
           fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
       } else if (lithGroupValue[j-1] == 1 &
           lithGroupValue[j + 1] == 1 &
           lithGroupThick[j - 1] > lithGroupThick[j + 1] ) {

           ecsDepth[k] = lithGroupDepth[j];
           fprintf(stderr, "shale break at  = %f \n ", ecsDepth[k]);
           k++;

       }
    }

   numECS = k + 1;

   // generate name of element sets
   grpcnt = 1;
   for (k = numECS - 2; k >= 0; k--, grpcnt++) {
       sprintf(ecsName[k], "ECS-%d", grpcnt);
       if (grpcnt % 2 == 0) {
           strncpy(ecsColour[k], "cyan", 4);
       }
       else {
           strncpy(ecsColour[k], "gray", 4);
       }
   }

   //write out results to .csv file
   FILE* ecsfile;
   ecsfile = fopen("./data/ecs.csv", "w");
   if (ecsfile != NULL) {
       fprintf(stderr, "./data/ecsfile.csv has been opened for write \n");
   }

   fprintf(stderr, "number of ECS = %d \n", numECS);

   for (k = 0; k < numECS; k++) {
       fprintf(ecsfile, "%f, %s, %s\n", ecsDepth[k], ecsName[k], ecsColour[k]);
   }

   fclose(ecsfile);

    


}
