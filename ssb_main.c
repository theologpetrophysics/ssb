

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

void calcBedVsh(
    int numLogPts,
    double* valueLogDepth,
    double* valueLogVsh,
    int numBeds,
    double* bedDepth,
    double* bedThick,
    double* bedVsh,
    double* bedTSF)
{

    int i, j, vshcnt;
    double vshsum;

    for (j = 0; j < numBeds; j++) {
        // calculate bed thickness
        bedThick[j] = bedDepth[j + 1] - bedDepth[j];

        // caluclate bed average VSH
        vshsum = 0;
        vshcnt = 0;
        for (i = 0; i < numLogPts; i++) {

            if (valueLogDepth[i] >= bedDepth[j] &
                valueLogDepth[i] < bedDepth[j + 1]) {

                vshsum += valueLogVsh[i];
                vshcnt++;
            }
        }
        bedVsh[j] = vshsum / (double)vshcnt;

        bedTSF[j] = bedThick[j] / (1.0 - bedVsh[j]);
    }

    // caluclate Vsh


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
    double lithGroupMinThick)

{
    int i, k, j, grpcnt;

    double smthData[1000] = { 0 }, localsValue;
    int framesToBeSmth;



    int debug = 0;

    // inform the user ******************
    fprintf(stderr, "Number of log samples to process is : %i \n The smoothfactor halfwindow is : %i \n", numpts, halfSmthWindow);
    //***********************************

    // Calculate smoothed VSH curve *****
    fprintf(stderr, "Generating smoothed Vsh curve... \n \n");
    //***********************************

    // logic to smooth the VSH curve
    //Cycle through data frames
    for (i = 0; i < numpts; i++) {

        if (halfSmthWindow == 0) {
            vshSmth[i] = vsh[i];

        }
        else {
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
    }

    // caluclate first and second derivative of Vsh and smoothed Vsh
    for (i = 1; i < numpts; i++) {
        vshFirstDeriv[i] = (vsh[i] - vsh[i - 1]) / (depth[i] - depth[i - 1]);
        vshSmthFirstDeriv[i] = (vshSmth[i] - vshSmth[i - 1]) / (depth[i] - depth[i - 1]);
    }

    for (i = 1; i < numpts; i++) {
        vshSecondDeriv[i] = (vshFirstDeriv[i] - vshFirstDeriv[i -1 ]) / (depth[i] - depth[i - 1]);
        vshSmthSecondDeriv[i] = (vshSmthFirstDeriv[i] - vshSmthFirstDeriv[i - 1]) / (depth[i] - depth[i - 1]);
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //--------------------------------------------------------------------------------------------
    // logic to calculate elements (bedsets)-------------------------------------

    int numE;
    int eframecnt;
    double ethicksum, ethickmin, ethickmax, evshsum;
    double* eGroupDepth;
    double* eGroupThick;
    double* eGroupVsh;
    double* eGroupValue;
    double* eGroupTSF;
    allocateMemory1DD(&eGroupDepth, numpts, 0);
    allocateMemory1DD(&eGroupThick, numpts, 0);
    allocateMemory1DD(&eGroupVsh, numpts, 0);
    allocateMemory1DD(&eGroupTSF, numpts, 0);
    allocateMemory1DD(&eGroupValue, numpts, 0);
    char eGroupName[LGROUPSIZE][10] = { {0} };
    char eGroupLith[LGROUPSIZE][5] = { {0} };
    char esGroupColour[LGROUPSIZE][15] = { {0} };


    // First pick sand/shale depths
    for (i = 0; i < numpts; i++) {

        if (vshSmth[i] == 1.0) {
            lithValue[i] = 1;
        }
        else if (vshSmth[i] == 0) {
            lithValue[i] = 0;
        }
        // use average vsh method for lithogroups
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
            // use derivative method for lithogroups
            if (vshSmthSecondDeriv[i] >= 0) {
                lithValue[i] = 0;
            }
            else {
                lithValue[i] = 1;
            }
        }
    }

    //lithgroup counter
    j = 0;

    // start with top frame and define group tops
    eGroupDepth[j] = depth[0];
    if (lithValue[0] == 1) {
        eGroupValue[j] = 1;
    }
    else {
        eGroupValue[j] = 0;
    }

    j++;

    for (i = 1; i < numpts-1; i++) {
        if (lithValue[i] != lithValue[i - 1] &
            lithValue[i] == lithValue[i+1] ) {

            // define current element
            eGroupValue[j] = lithValue[i];
            eGroupDepth[j] = depth[i];
            
            j++;

            if (debug == 1) {
                fprintf(stderr, "grp depth  = %f \n ", eGroupDepth[j]);
            }

        }
    }

    i = numpts-1;
    // last frame of log data
    eGroupValue[j] = MISSING;
    eGroupDepth[j] = depth[i];
    eGroupThick[j] = MISSING;
    eGroupVsh[j] = MISSING;
    j++;
   
    numE = j - 1;

    // generate thickness, VSH and TSF of each ES
    calcBedVsh(
        numpts,
        depth,
        vsh,
        numE,
        eGroupDepth,
        eGroupThick,
        eGroupVsh,
        eGroupTSF);

/*
    fprintf(stderr, "Number of Element Sets = %d \n", numE);
    fprintf(stderr, "Average thickness of Element Sets = %f \n", ethicksum / numE);
    fprintf(stderr, "Minimum thickness of Element Sets = %f \n", ethickmin);
    fprintf(stderr, "Maximum thickness of Element Sets = %f \n", ethickmax);
    */

    // generate lithology of Element sets
    for (i = 0; i < numE; i++) {
        if (eGroupValue[i] == 1) {
            //fprintf(stderr, "lithgroup %d - ", i);
            strncpy(eGroupLith[i], "sh", 2);
            strncpy(esGroupColour[i], "gray15", 6);
        }
        else {
            strncpy(eGroupLith[i], "ss", 2);
            strncpy(esGroupColour[i], "yellow", 6);
        }
    }
       
    // generate name of element sets
    grpcnt = 1;
    for (i = numE - 1; i >= 0; i--, grpcnt++) {
        sprintf(eGroupName[i], "E-%d", grpcnt);        
    }
  
    //write out results to .csv file
    FILE* vshlgfile;

    //fprintf(stderr, "numgrps = %d \n", numE);

    vshlgfile = fopen("./data/elements.csv", "w");
    if (vshlgfile != NULL) {
        fprintf(stderr, "./data/elements.csv has been opened for write \n");
    }

    fprintf(stderr, "numer of Elements = %d \n", numE);

    grpcnt = 1;
    for (i = numE - 1; i >= 0; i--) {
        //fprintf(stderr, "name = %s \n", lithGroupName[i]);
    }

    for (i = 0; i <= numE; i++) {
        fprintf(vshlgfile, "%f, %f, %s, %s, %s, %f \n", eGroupDepth[i], eGroupThick[i], eGroupLith[i], eGroupName[i], esGroupColour[i], eGroupVsh[i]);
    }

    fclose(vshlgfile);

    //  end of potential sub function for generating the element set data--------------------------------------
    //---------------------------------------------------------------------------------------------------------


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // logic to calculate Element Sets (ES)---------------------
    //----------------------------------------------------------

    int numES;
    int esframecnt;
    double esthicksum, esthickmin, esthickmax, esvshsum;
    double* esDepth;
    double* esVsh;
    double* esThick;
    double* esTSF;
    allocateMemory1DD(&esDepth, numE, 0);
    allocateMemory1DD(&esVsh, numE, 0);
    allocateMemory1DD(&esThick, numE, 0);
    allocateMemory1DD(&esTSF, numE, 0);
    char esName[LGROUPSIZE][10] = { {0} };
    char esColour[LGROUPSIZE][15] = { {0} };

    j = 0;

    for (i = 0; i <= numE; i++) {
        if (i == 0) {
            // top surface is treated as top of ES
            esDepth[j] = eGroupDepth[i];
            j++;
        }
        else if (i == numE) {
            // base surface is base of ES
            esDepth[j] = eGroupDepth[numE];
            j++;
         }
        else if (eGroupValue[i] == 0 &
            eGroupValue[i-1] == 1 &
            eGroupValue[i+1] == 1) {
            // top of sand/shale couplet is top of ES unit
            esDepth[j] = eGroupDepth[i];
            j++;
        }
    }

    numES = j - 1;

    // generate name of Element Sets
    grpcnt = 1;
    for (k = numES - 1; k >= 0; k--, grpcnt++) {
        sprintf(esName[k], "ES-%d", grpcnt);
        if (grpcnt % 2 == 0) {
            strncpy(esColour[k], "light_cyan", 10);
        }
        else {
            strncpy(esColour[k], "light_gray", 10);
        }
    }

    // generate thickness, VSH and TSF of each ES
    calcBedVsh(
        numpts,
        depth,
        vsh,
        numES,
        esDepth,
        esThick,
        esVsh,
        esTSF);


    //write out results to .csv file
    FILE* esfile;
    esfile = fopen("./data/es.csv", "w");
    if (esfile != NULL) {
        fprintf(stderr, "./data/esfile.csv has been opened for write \n");
    }

    fprintf(stderr, "number of Element Sets = %d \n", numES);

    for (k = 0; k <= numES; k++) {
        fprintf(esfile, "%f, %s, %s, %f, %f, %f \n", esDepth[k], esName[k], esColour[k], esThick[k], esVsh[k], esTSF[k]);
    }

    fclose(esfile);
    

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  Now need to find the shale breaks....
    //  Find when shale is thicker above than below
    //  if shale above current surface is thicker than next shale below then current surface is a shale break
    //  Zone between shale breaks is called Element Complex Set (ECS)

    // top surface is treated as top of ECS

    int numECS;
    int ecsframecnt;
    double ecsthicksum, ecsthickmin, ecsthickmax, ecsvshsum;
    double* ecsDepth;
    double* ecsVsh;
    double* ecsThick;
    double* ecsTSF;
    allocateMemory1DD(&ecsDepth, numES, 0);
    allocateMemory1DD(&ecsVsh, numES, 0);
    allocateMemory1DD(&ecsThick, numES, 0);
    allocateMemory1DD(&ecsTSF, numES, 0);
    char ecsName[LGROUPSIZE][10] = { {0} };
    char ecsColour[LGROUPSIZE][15] = { {0} };

    k = 0;
   /* ecsthicksum = 0;
    ecsthickmin = 100;
    ecsthickmax = 0;
    ecsvshsum = 0;
  */


   for (j = 0; j <= numE; j++) {
       if (j == 0) {
           ecsDepth[k] = eGroupDepth[j];
           k++;
           //fprintf(stderr, "Top at  = %f \n ", ecsDepth[k]);
           
       }
       else if (j == numE) {
           ecsDepth[k] = eGroupDepth[j];
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       } else if (eGroupValue[j-1] == 1 &
           eGroupValue[j + 1] == 1 &
           eGroupThick[j - 1] > eGroupThick[j + 1] ) {
           ecsDepth[k] = eGroupDepth[j];
           k++;
           //fprintf(stderr, "shale break at  = %f \n ", ecsDepth[k]);
           
       }
    }

   numECS = k - 1;

   // generate name of Element Complex Sets
   grpcnt = 1;
   for (k = numECS - 1; k >= 0; k--, grpcnt++) {
       sprintf(ecsName[k], "ECS-%d", grpcnt);
       if (grpcnt % 2 == 0) {
           strncpy(ecsColour[k], "tan", 3);
       }
       else {
           strncpy(ecsColour[k], "thistle", 7);
       }
   }



   // generate thickness, VSH and TSF of each ECS
   calcBedVsh(
       numpts,
       depth,
       vsh,
       numECS,
       ecsDepth,
       ecsThick,
       ecsVsh,
       ecsTSF);
       

   //write out results to .csv file
   FILE* ecsfile;
   ecsfile = fopen("./data/ecs.csv", "w");
   if (ecsfile != NULL) {
       fprintf(stderr, "./data/ecsfile.csv has been opened for write \n");
   }

   fprintf(stderr, "number of Element Complex Sets = %d \n", numECS);

   for (k = 0; k <= numECS; k++) {
       fprintf(ecsfile, "%f, %s, %s, %f, %f, %f \n", ecsDepth[k], ecsName[k], ecsColour[k], ecsThick[k], ecsVsh[k], ecsTSF[k]);
   }

   fclose(ecsfile);



   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  Now need to find parasequence tops....
   //  Find when TSF goes up
   //  Zone between increasinf TSF tops are called parasequences

   int numPS;
   int psframecnt;
   double psthicksum, psthickmin, psthickmax, psvshsum;
   double* psDepth;
   double* psVsh;
   double* psThick;
   double* psTSF;
   allocateMemory1DD(&psDepth, numECS, 0);
   allocateMemory1DD(&psVsh, numECS, 0);
   allocateMemory1DD(&psThick, numECS, 0);
   allocateMemory1DD(&psTSF, numECS, 0);
   char psName[LGROUPSIZE][10] = { {0} };
   char psColour[LGROUPSIZE][20] = { {0} };

   k = 0;

   for (j = 0; j <= numECS; j++) {
       if (j == 0) {
           psDepth[k] = ecsDepth[j];
           k++;
           //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

       }
       else if (j == numECS) {
           psDepth[k] = ecsDepth[j];
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       }
       else if (ecsTSF[j-1] > ecsTSF[j]) {
           psDepth[k] = ecsDepth[j];
           k++;
       }
   }

   numPS = k - 1;

   // generate name of Element Complex Sets
   grpcnt = 1;
   for (k = numPS - 1; k >= 0; k--, grpcnt++) {
       sprintf(psName[k], "PS-%d", grpcnt);
       if (grpcnt % 2 == 0) {
           strncpy(psColour[k], "yellow_green", 12);
       }
       else {
           strncpy(psColour[k], "light_sea_green", 15);
       }
   }

   // generate thickness, VSH and TSF of each ECS
   calcBedVsh(
       numpts,
       depth,
       vsh,
       numPS,
       psDepth,
       psThick,
       psVsh,
       psTSF);


   //write out results to .csv file
   FILE* psfile;
   psfile = fopen("./data/ps.csv", "w");
   if (psfile != NULL) {
       fprintf(stderr, "./data/ps.csv has been opened for write \n");
   }

   fprintf(stderr, "number of Para-sequences = %d \n", numPS);

   for (k = 0; k <= numPS; k++) {
       fprintf(psfile, "%f, %s, %s, %f, %f, %f \n", psDepth[k], psName[k], psColour[k], psThick[k], psVsh[k], psTSF[k]);
   }

   fclose(psfile);





    


}
