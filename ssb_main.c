

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "theolog_general_functions.h"

#define MISSING -999.25
#define LGROUPSIZE 2000


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

        bedTSF[j] = bedThick[j] / (1.0 - fmin(bedVsh[j],0.95));
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
    char* ecsMethod,
    double ecsShaleBreakLimit,
    char* tsfPickMethod,
    char* optElementLog,
    double* vshSmth,
    double* grSmth,
    double* firstDeriv,
    double* secondDeriv,
    double* lithValue,
    double lithGroupMinThick)

{
    int i, k, j, m, grpcnt;

    double smthData[1000] = { 0 }, localsValue;
    double grsmthData[1000] = { 0 };
    int framesToBeSmth;



    int debug = 0;

    // inform the user ******************
    fprintf(stderr, "Number of log samples to process is : %i \n The smoothfactor halfwindow is : %i \n", numpts, halfSmthWindow);
    //***********************************
    fprintf(stderr, "before deriv... \n \n");
    // Calculate smoothed VSH curve *****
    fprintf(stderr, "Generating smoothed Vsh curve... \n \n");
    //***********************************

    // logic to smooth the VSH curve
    //Cycle through data frames
    for (i = 0; i < numpts; i++) {

        if (halfSmthWindow == 0) {
            vshSmth[i] = vsh[i];
            grSmth[i] = gr[i];

        }
        else {
            j = 0;
            //generate array of data to be averaged
            for (k = fmax(i - halfSmthWindow, 0); k <= fmin((i + halfSmthWindow), numpts - 1); k++, j++) {
                smthData[j] = vsh[k];
                grsmthData[j] = gr[k];
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
            localsValue = smoothLogData(framesToBeSmth, "box", grsmthData);
            grSmth[i] = localsValue;
            if (debug == 1) {
                fprintf(stderr, "Smoothed value = %f \n ", localsValue);
            }
        }
    }

    fprintf(stderr, "before deriv... \n \n");

    // use GR or VSH for the derivative method
    double* derivInLog;
    allocateMemory1DD(&derivInLog, numpts, 0);

    if (strncmp(optElementLog, "GR", 2) == 0) {
        for (i = 0; i < numpts; i++) {
            derivInLog[i] = grSmth[i];
        }
    }
    else if (strncmp(optElementLog, "VSH", 3) == 0) {
        for (i = 0; i < numpts; i++) {
            derivInLog[i] = vshSmth[i];
        }
    }

    // calculate first and second derivative of input log
    for (i = 0; i < numpts; i++) {
        firstDeriv[i] = (derivInLog[i+1] - derivInLog[i]) / (depth[i+1] - depth[i]);
    }

    for (i = 1; i < numpts; i++) {
        secondDeriv[i] = (firstDeriv[i + 1] - firstDeriv[i]) / (depth[i + 1] - depth[i]);
    }

    fprintf(stderr, "before bedsets... \n \n");

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //--------------------------------------------------------------------------------------------
    // logic to calculate elements (bedsets)-------------------------------------

    int numE;
    int eframecnt;
    double ethicksum, ethickmin, ethickmax, evshsum;
    double* eDepth;
    double* eThick;
    double* eVsh;
    double* eValue;
    double* eTSF;
    allocateMemory1DD(&eDepth, numpts, 0);
    allocateMemory1DD(&eThick, numpts, 0);
    allocateMemory1DD(&eVsh, numpts, 0);
    allocateMemory1DD(&eTSF, numpts, 0);
    allocateMemory1DD(&eValue, numpts, 0);
    char eName[LGROUPSIZE][10] = { {0} };
    char eLith[LGROUPSIZE][5] = { {0} };
    char eColour[LGROUPSIZE][15] = { {0} };


    // First pick sand/shale depths
    for (i = 0; i < numpts; i++) {

        if (vshSmth[i] == 1.0 & strncmp(optElementLog, "VSH", 3) == 0) {
            lithValue[i] = 1;
        }
        else if (vshSmth[i] == 0 & strncmp(optElementLog, "VSH", 3) == 0) {
            lithValue[i] = 0;
        }
        // use average vsh method for lithogroups
        else if (strncmp(lithGroupMethod, "AVERAGE", 7) == 0) {

            if (strncmp(optElementLog, "VSH", 3) == 0 &
                vsh[i] >= vshSmth[i] |
                strncmp(optElementLog, "GR", 2) == 0 &
                gr[i] >= grSmth[i]) {
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
            if (secondDeriv[i] >= 0) {
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
    eDepth[j] = depth[0];
    if (lithValue[0] == 1) {
        eValue[j] = 1;
    }
    else {
        eValue[j] = 0;
    }

    j++;

    for (i = 1; i < numpts-1; i++) {
        /*if (lithValue[i] != lithValue[i - 1] &
            lithValue[i] == lithValue[i+1] ) { */
        if (lithValue[i] != lithValue[i - 1]) {

            // define current element
            eValue[j] = lithValue[i];
            eDepth[j] = depth[i];
            
            j++;

            if (debug == 1) {
                fprintf(stderr, "grp depth  = %f \n ", eDepth[j]);
            }

        }
    }

    i = numpts-1;
    // last frame of log data
    eValue[j] = MISSING;
    eDepth[j] = depth[i];
    eThick[j] = MISSING;
    eVsh[j] = MISSING;
    j++;
   
    numE = j - 1;

    // generate thickness, VSH and TSF of each ES
    calcBedVsh(
        numpts,
        depth,
        vsh,
        numE,
        eDepth,
        eThick,
        eVsh,
        eTSF);


    // generate lithology of Element sets
    for (i = 0; i < numE; i++) {
        if (eValue[i] == 1) {
            //fprintf(stderr, "lithgroup %d - ", i);
            strncpy(eLith[i], "sh", 2);
            strncpy(eColour[i], "gray", 6);
        }
        else {
            strncpy(eLith[i], "ss", 2);
            strncpy(eColour[i], "yellow", 6);
        }
    }
       
    // generate name of element sets
    grpcnt = 1;
    for (i = numE - 1; i >= 0; i--, grpcnt++) {
        sprintf(eName[i], "E-%d", grpcnt);        
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
        fprintf(vshlgfile, "%f, %f, %s, %s, %s, %f \n", eDepth[i], eThick[i], eLith[i], eName[i], eColour[i], eVsh[i]);
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
    double* esShaleThick;
    allocateMemory1DD(&esDepth, numE, 0);
    allocateMemory1DD(&esVsh, numE, 0);
    allocateMemory1DD(&esThick, numE, 0);
    allocateMemory1DD(&esTSF, numE, 0);
    allocateMemory1DD(&esShaleThick, numE, 0);
    char esName[LGROUPSIZE][10] = { {0} };
    char esColour[LGROUPSIZE][15] = { {0} };

    j = 0;

    for (i = 0; i <= numE; i++) {
        if (i == 0) {
            // top surface is treated as top of ES
            esDepth[j] = eDepth[i];
            j++;
        }
        else if (i == numE) {
            // base surface is base of ES
            esDepth[j] = eDepth[numE];
            j++;
         }

        else if (i == numE-1 & eValue[i] == 0) {
            // base surface is base of ES
            esDepth[j] = eDepth[numE];
            j++;
        }
        else if (eValue[i] == 0 &
            eValue[i-1] == 1 &
            eValue[i+1] == 1) {
            // top of sand/shale couplet is top of ES unit
            esDepth[j] = eDepth[i];
            esShaleThick[j] = eThick[i + 1];
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
            strncpy(esColour[k], "light_blue", 10);
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
        fprintf(esfile, "%f, %s, %s, %f, %f, %f, %f \n", esDepth[k], esName[k], esColour[k], esThick[k], esShaleThick[k], esVsh[k], esTSF[k]);
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
    double* ecsSbThick;
    double* ecsSsThick;
    double* ecsTSF;
    double* ecsESsVsh;
    double* ecsEShVsh;
    double* ecsSbRatio;

    allocateMemory1DD(&ecsDepth, numES, 0);
    allocateMemory1DD(&ecsVsh, numES, 0);
    allocateMemory1DD(&ecsThick, numES, 0);
    allocateMemory1DD(&ecsSbThick, numES, 0);
    allocateMemory1DD(&ecsSbRatio, numES, 0);
    allocateMemory1DD(&ecsSsThick, numES, 0);
    allocateMemory1DD(&ecsTSF, numES, 0);
    allocateMemory1DD(&ecsESsVsh, numES, 0);
    allocateMemory1DD(&ecsEShVsh, numES, 0);
    char ecsName[LGROUPSIZE][10] = { {0} };
    char ecsColour[LGROUPSIZE][15] = { {0} };

    k = 0;
   /* ecsthicksum = 0;
    ecsthickmin = 100;
    ecsthickmax = 0;
    ecsvshsum = 0;
  */

    if (strncmp(ecsMethod, "SB", 2) == 0) {

        for (j = 0; j <= numES; j++) {
            if (j == 0) {
                ecsDepth[k] = esDepth[j];
                k++;
                //fprintf(stderr, "Top at  = %f \n ", ecsDepth[k]);

            }
            else if (j == numES) {
                ecsDepth[k] = esDepth[j];
                //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
                k++;

            }
            else if (esShaleThick[j - 1] > esShaleThick[j] * (1.0 + ecsShaleBreakLimit)) {
                
                ecsDepth[k] = esDepth[j];
                ecsSbThick[k] = esShaleThick[j - 1];
                ecsSsThick[k] = esThick[j - 1] - esShaleThick[j - 1];
                // calculate shale break thickness change ratio for calibration
                ecsSbRatio[k] = esShaleThick[j - 1] / esShaleThick[j];
                k++;
                
                //fprintf(stderr, "shale thick (%s)  = %f \n ", esName[j - 1], esShaleThick[j - 1]);
                //fprintf(stderr, "shale thick (%s)  = %f \n ", esName[j], esShaleThick[j]);
                //fprintf(stderr, "---\n"); 
            }
        }
    }
    else if (strncmp(ecsMethod, "TSF", 3) == 0) {
        for (j = 0; j <= numES; j++) {
            if (j == 0) {
                ecsDepth[k] = esDepth[j];
                k++;
                //fprintf(stderr, "Top at  = %f \n ", ecsDepth[k]);

            }
            else if (j == numES) {
                ecsDepth[k] = esDepth[j];
                //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
                k++;

            }
            else if (esTSF[j] < esTSF[j + 1] &
                esTSF[j - 1] > esTSF[j] &
                esTSF[j - 2] < esTSF[j - 1]) {

                ecsDepth[k] = esDepth[j];
                fprintf(stderr, "ecs type 1  = %f \n ", ecsDepth[k]);
                k++;

            } else if (esTSF[j + 1] < esTSF[j] &
                esTSF[j + 2] > esTSF[j + 1] &
                esTSF[j - 1] > esTSF[j] &
                esTSF[j - 2] < esTSF[j - 1]) {

                ecsDepth[k] = esDepth[j];
                fprintf(stderr, "ecs type 2  = %f \n ", ecsDepth[k]);
                k++;

            } else if (esTSF[j - 1] > esTSF[j] &
                esTSF[j - 2] < esTSF[j - 1]) {

            ecsDepth[k] = esDepth[j];
            fprintf(stderr, "ecs type 3  = %f \n ", ecsDepth[k]);
            k++;

            }
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

   // calculate the average Vsh of SS/SH elements within each ECS
   double ssvshsum, ssthicksum, shvshsum, shthicksum;

   for (j = 0; j <= numECS; j++) {
       ssvshsum = 0;
       ssthicksum = 0;
       shvshsum = 0;
       shthicksum = 0;

       for (i = 0; i <= numE; i++) {
           if (eDepth[i] >= ecsDepth[j] &
               eDepth[i] < ecsDepth[j + 1])
           {
               if (eValue[i] == 0)
               {
                   ssvshsum += eVsh[i] * eThick[i];
                   ssthicksum += eThick[i];
               }
               ecsESsVsh[j] = ssvshsum / ssthicksum;

               if (eValue[i] == 1)
               {
                   shvshsum += eVsh[i] * eThick[i];
                   shthicksum += eThick[i];
               }
               ecsEShVsh[j] = shvshsum / shthicksum;
           }
       }
   }
       

   //write out results to .csv file
   FILE* ecsfile;
   ecsfile = fopen("./data/ecs.csv", "w");
   if (ecsfile != NULL) {
       fprintf(stderr, "./data/ecsfile.csv has been opened for write \n");
   }

   fprintf(stderr, "number of Element Complex Sets = %d \n", numECS);

   for (k = 0; k <= numECS; k++) {
       fprintf(ecsfile, "%f, %s, %s, %f, %f, %f, %f, %f, %f, %f, %f \n", ecsDepth[k], ecsName[k], ecsColour[k], 
           ecsThick[k], ecsVsh[k], ecsTSF[k], ecsSbThick[k], ecsSbRatio[k], ecsSsThick[k], ecsESsVsh[k], ecsEShVsh[k]);
   }

   fclose(ecsfile);

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  Now need to find parasequence tops....
   //  Find when TSF goes up (and VSH goes up)
   //  Zone between increasinf TSF tops are called parasequences

   
   //---------------------------------------------------------------------------------
   // find basic PS tops
   int numBasicPS;
   double* psBasicDepth;
   allocateMemory1DD(&psBasicDepth, numECS, 0);
   char psBasicMethod[LGROUPSIZE][20] = { {0} };
   char psBasicName[LGROUPSIZE][10] = { {0} };
   k = 0;
   m = 0;

   for (j = 0; j <= numECS; j++) {
       if (j == 0) {
           psBasicDepth[k] = ecsDepth[j];
           strncpy(psBasicMethod[k], "Top", 15);
           k++;
           //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

       }
       else if (j == numECS) {
           psBasicDepth[k] = ecsDepth[j];
           strncpy(psBasicMethod[k], "Base", 15);
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       }
       // very basic TSF increase (local TSF minimum)
       else if (ecsTSF[j - 1] > ecsTSF[j] * 1.05)
       {
           psBasicDepth[k] = ecsDepth[j];
           strncpy(psBasicMethod[k], "TSf only", 15);
           k++;
       }
   }

   numBasicPS = k - 1;

   // generate name of Element Complex Sets
   grpcnt = 1;
   for (k = numBasicPS - 1; k >= 0; k--, grpcnt++) {
       sprintf(psBasicName[k], "PS-%d", grpcnt);
   }

   //write out results to .csv file
   FILE* psBasicfile;
   psBasicfile = fopen("./data/ps_basic.csv", "w");
   if (psBasicfile != NULL) {
       fprintf(stderr, "./data/ps_basic.csv has been opened for write \n");
   }

   fprintf(stderr, "number of Para-sequences = %d \n", numBasicPS);

   for (k = 0; k <= numBasicPS; k++) {
       fprintf(psBasicfile, "%f, %s, %s \n", psBasicDepth[k], psBasicName[k], psBasicMethod[k]);
   }

   fclose(psBasicfile);


   //---------------------------------------------------------------------------------
   // find basic PS tops
   int numTsfVshPS;
   double* psTsfVshDepth;
   allocateMemory1DD(&psTsfVshDepth, numECS, 0);
   char psTsfVshMethod[LGROUPSIZE][20] = { {0} };
   char psTsfVshName[LGROUPSIZE][10] = { {0} };
   k = 0;

   for (j = 0; j <= numECS; j++) {
       if (j == 0) {
           psTsfVshDepth[k] = ecsDepth[j];
           strncpy(psTsfVshMethod[k], "Top", 15);
           k++;
           //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

       }
       else if (j == numECS) {
           psTsfVshDepth[k] = ecsDepth[j];
           strncpy(psTsfVshMethod[k], "Base", 15);
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       }
       // very basic TSF increase (local TSF minimum)
       else if (ecsTSF[j - 1] > ecsTSF[j] * 1.05 &
                ecsVsh[j - 1] >= ecsVsh[j])
       {
           psTsfVshDepth[k] = ecsDepth[j];
           strncpy(psTsfVshMethod[k], "TSf/Vsh", 15);
           k++;
       }
   }

   numTsfVshPS = k - 1;

   // generate name of Element Complex Sets
   grpcnt = 1;
   for (k = numTsfVshPS - 1; k >= 0; k--, grpcnt++) {
       sprintf(psTsfVshName[k], "PS-%d", grpcnt);
   }

   //write out results to .csv file
   FILE* psTsfVshfile;
   psTsfVshfile = fopen("./data/ps_tsfvsh.csv", "w");
   if (psTsfVshfile != NULL) {
       fprintf(stderr, "./data/ps_tsfvsh.csv has been opened for write \n");
   }

   fprintf(stderr, "number of Para-sequences = %d \n", numTsfVshPS);

   for (k = 0; k <= numTsfVshPS; k++) {
       fprintf(psTsfVshfile, "%f, %s, %s \n", psTsfVshDepth[k], psTsfVshName[k], psTsfVshMethod[k]);
   }

   fclose(psTsfVshfile);
   
   //---------------------------------------------------------------------------------
   // find tsf trend PS tops
   int numTrendPS;
   double* psTrendDepth;
   allocateMemory1DD(&psTrendDepth, numECS, 0);
   char psTrendMethod[LGROUPSIZE][30] = { {0} };
   char psTrendName[LGROUPSIZE][20] = { {0} };
   k = 0;
   
   for (j = 0; j <= numECS; j++) {
       if (j == 0) {
           psTrendDepth[k] = ecsDepth[j];
           strncpy(psTrendMethod[k], "Top", 15);
           k++;
           //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

       }
       else if (j == numECS) {
           psTrendDepth[k] = ecsDepth[j];
           strncpy(psTrendMethod[k], "base", 4);
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       }
       
       // TSF increase with decreasing trend above
       else if (ecsTSF[j - 1] > ecsTSF[j] &
           ecsVsh[j - 1] >= ecsVsh[j] &
           ecsTSF[j - 2] < ecsTSF[j - 1])
       {
           psTrendDepth[k] = ecsDepth[j];
           strncpy(psTrendMethod[k], "TSF/Vsh + up trend", 19);
           k++;
       }
   }
   
   
   numTrendPS = k - 1;
   fprintf(stderr, "number of Trend Para-sequences = %d \n", numTrendPS);
   
   // generate name of Element Complex Sets
   grpcnt = 1;
   
   for (k = numTrendPS - 1; k >= 0; k--, grpcnt++) {
       sprintf(psTrendName[k], "PS-%d", grpcnt);
   }
   
  
   //write out results to .csv file
   FILE* psTrendfile;
   psTrendfile = fopen("./data/ps_trend.csv", "w");
   if (psTrendfile != NULL) {
       fprintf(stderr, "./data/ps_trend.csv has been opened for write \n");
   }


   for (k = 0; k <= numTrendPS; k++) {
       fprintf(psTrendfile, "%f, %s, %s \n", psTrendDepth[k], psTrendName[k], psTrendMethod[k]);
   }

   fclose(psTrendfile);

   
   
   //---------------------------------------------------------------------------------
  // find tsf transgressive PS tops
   int numTransPS;
   double* psTransDepth;
   allocateMemory1DD(&psTransDepth, numECS, 0);
   char psTransMethod[LGROUPSIZE][20] = { {0} };
   char psTransName[LGROUPSIZE][10] = { {0} };
   k = 0;

   for (j = 0; j <= numECS; j++) {
       if (j == 0) {
           psTransDepth[k] = ecsDepth[j];
           strncpy(psTransMethod[k], "Top", 15);
           k++;
           //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

       }
       else if (j == numECS) {
           psTransDepth[k] = ecsDepth[j];
           strncpy(psTransMethod[k], "Base", 15);
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       }
       // allowance for transgressive top
       else if ((strncmp(tsfPickMethod, "TSF", 3) == 0 |
           strncmp(tsfPickMethod, "EITHER", 6) == 0) &
           ecsTSF[j - 1] > ecsTSF[j] &
           ecsVsh[j - 1] >= ecsVsh[j] &
           ecsTSF[j - 2] < ecsTSF[j - 1] &
           ecsTSF[j + 1] < ecsTSF[j] &
           ecsTSF[j + 2] > ecsTSF[j + 1])

       {
           psTransDepth[k] = ecsDepth[j];
           strncpy(psTransMethod[k], "TSF/Vsh with ts ", 20);
           k++;
       }
   }

   numTransPS = k - 1;

   // generate name of Element Complex Sets
   grpcnt = 1;
   for (k = numTransPS - 1; k >= 0; k--, grpcnt++) {
       sprintf(psTransName[k], "PS-%d", grpcnt);
   }

   //write out results to .csv file
   FILE* psTransfile;
   psTransfile = fopen("./data/ps_trans.csv", "w");
   if (psTransfile != NULL) {
       fprintf(stderr, "./data/ps_trans.csv has been opened for write \n");
   }

   fprintf(stderr, "number of Para-sequences = %d \n", numTransPS);

   for (k = 0; k <= numTransPS; k++) {
       fprintf(psTransfile, "%f, %s, %s \n", psTransDepth[k], psTransName[k], psTransMethod[k]);
   }

   fclose(psTransfile);

   //---------------------------------------------------------------------------------
  // find SB thickness tops
   int numSbPS;
   double* psSbDepth;
   allocateMemory1DD(&psSbDepth, numECS, 0);
   char psSbMethod[LGROUPSIZE][20] = { {0} };
   char psSbName[LGROUPSIZE][10] = { {0} };
   k = 0;

   for (j = 0; j <= numECS; j++) {
       if (j == 0) {
           psSbDepth[k] = ecsDepth[j];
           strncpy(psSbMethod[k], "Top", 15);
           k++;
           //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

       }
       else if (j == numECS) {
           psSbDepth[k] = ecsDepth[j];
           strncpy(psSbMethod[k], "Base", 15);
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       }
       else if (ecsSbThick[j] > ecsSbThick[j + 1] * 1.25 &
                ecsESsVsh[j] < ecsESsVsh[j - 1])
       {
           psSbDepth[k] = ecsDepth[j];
           strncpy(psSbMethod[k], "SB thickness & GS", 17);
           k++;
       }
   }

   numSbPS = k - 1;

   // generate name of Element Complex Sets
   grpcnt = 1;
   for (k = numSbPS - 1; k >= 0; k--, grpcnt++) {
       sprintf(psSbName[k], "PS-%d", grpcnt);
   }

   //write out results to .csv file
   FILE* psSbfile;
   psSbfile = fopen("./data/ps_sb.csv", "w");
   if (psSbfile != NULL) {
       fprintf(stderr, "./data/ps_sb.csv has been opened for write \n");
   }

   fprintf(stderr, "number of Para-sequences = %d \n", numSbPS);

   for (k = 0; k <= numSbPS; k++) {
       fprintf(psSbfile, "%f, %s, %s \n", psSbDepth[k], psSbName[k], psSbMethod[k]);
   }

   fclose(psSbfile);
   
   //-----------------------------------------------------------------------------------
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
   char psMethod[LGROUPSIZE][30] = { {0} };

   k = 0;
   for (j = 0; j <= numECS; j++) {
       if (j == 0) {
           psDepth[k] = ecsDepth[j];
           strncpy(psMethod[k], "Top", 15);
           k++;
           //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

       }
       else if (j == numECS) {
           psDepth[k] = ecsDepth[j];
           strncpy(psMethod[k], "Base", 15);
           //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
           k++;

       }
       // PS picking criteria
       else if ((strncmp(tsfPickMethod, "TSF", 3) == 0 |
           strncmp(tsfPickMethod, "EITHER", 6) == 0) &
           ecsTSF[j - 1] > ecsTSF[j] * 1.05)
       {

           if (ecsVsh[j - 1] >= ecsVsh[j] &
               ecsTSF[j - 2] < ecsTSF[j - 1])
           {
               psDepth[k] = ecsDepth[j];
               strncpy(psMethod[k], "TSF-VSH-TREND", 25);
               k++;
           }

           else if (ecsVsh[j - 1] >= ecsVsh[j])
           {
               psDepth[k] = ecsDepth[j];
               strncpy(psMethod[k], "TSF-VSH", 15);
               k++;
           }
           else if (ecsTSF[j - 2] < ecsTSF[j - 1])
           {
               psDepth[k] = ecsDepth[j];
               strncpy(psMethod[k], "TSF-TREND ", 20);
               k++;
           }
       }
       else if ((strncmp(tsfPickMethod, "SB", 2) == 0 |
               strncmp(tsfPickMethod, "EITHER", 6) == 0) &
               ecsSbThick[j] > ecsSbThick[j + 1] * 1.25 &
               ecsESsVsh[j] < ecsESsVsh[j - 1])
       {
           if (ecsVsh[j - 1] >= ecsVsh[j] &
               ecsTSF[j - 2] < ecsTSF[j - 1])
           {
               psDepth[k] = ecsDepth[j];
               strncpy(psMethod[k], "ShBr-VSH-TREND", 15);
               k++;
           }
           else if (ecsVsh[j - 1] >= ecsVsh[j])
           {
               psDepth[k] = ecsDepth[j];
               strncpy(psMethod[k], "ShBr-VSH", 15);
               k++;
           }
           else if (ecsTSF[j - 2] < ecsTSF[j - 1])
           {
               psDepth[k] = ecsDepth[j];
               strncpy(psMethod[k], "ShBr-TREND ", 20);
               k++;
           }
       }
       else if (ecsTSF[j - 1] > ecsTSF[j] * 1.05 &
                ecsSbThick[j] > ecsSbThick[j + 1] * 1.25 &
                ecsESsVsh[j] < ecsESsVsh[j - 1])
       {
            psDepth[k] = ecsDepth[j];
            strncpy(psMethod[k], "TSF-ShBr", 15);
            k++;
       }

       
       /*
       // TSF increase + VSH (local TSF minimum)
       else if ((strncmp(tsfPickMethod, "TSF", 3) == 0 |
           strncmp(tsfPickMethod, "EITHER", 6) == 0)  &
           ecsTSF[j - 1] > ecsTSF[j] * 1.05 &
           ecsVsh[j - 1] >= ecsVsh[j])
       {
           psDepth[k] = ecsDepth[j];
           strncpy(psMethod[k], "TSf/Vsh", 15);
           k++;
       } 
       // TSF increase with decreasing trend above
       else if ((strncmp(tsfPickMethod, "TSF", 3) == 0 | 
           strncmp(tsfPickMethod, "EITHER", 6) == 0) &
           ecsTSF[j - 1] > ecsTSF[j] &
           ecsVsh[j - 1] >= ecsVsh[j] &
           ecsTSF[j - 2] < ecsTSF[j - 1]) 
       {
           psDepth[k] = ecsDepth[j];
           strncpy(psMethod[k], "TSF/Vsh + up trend ", 18);
           k++;
       } 
        // allowance for transgressive top
       else if ((strncmp(tsfPickMethod, "TSF", 3) == 0 |
           strncmp(tsfPickMethod, "EITHER", 6) == 0) &
            ecsTSF[j - 1] > ecsTSF[j] &
            ecsVsh[j - 1] >= ecsVsh[j] &
            ecsTSF[j - 2] < ecsTSF[j - 1] &              
            ecsTSF[j + 1] < ecsTSF[j] &
            ecsTSF[j + 2] > ecsTSF[j + 1])
            
       {
           psDepth[k] = ecsDepth[j];
           strncpy(psMethod[k], "TSF/Vsh with ts ", 20);
           k++;
       }
       */       

       /*
       else if (strncmp(tsfPickMethod, "SB", 2) == 0 &
           ecsSbThick[j] >= ecsSbThick[j + 1] * 0.9 &
           ecsSsThick[j] < ecsSsThick[j + 1])
       {
           psDepth[k] = ecsDepth[j];
           strncpy(psMethod[k], "SS thickness", 12);
           k++;
       } */

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
       fprintf(psfile, "%f, %s, %s, %f, %f, %f, %s \n", psDepth[k], psName[k], psColour[k], psThick[k], psVsh[k], psTSF[k], psMethod[k]);
   }

   fclose(psfile);



   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  Now need to find fourth order sequences tops....
   //  

   int numFOS;
   double* fosDepth;
   double* fosVsh;
   double* fosThick;
   double* fosTSF;
   allocateMemory1DD(&fosDepth, numPS, 0);
   allocateMemory1DD(&fosVsh, numPS, 0);
   allocateMemory1DD(&fosThick, numPS, 0);
   allocateMemory1DD(&fosTSF, numPS, 0);
   char fosName[LGROUPSIZE][10] = { {0} };
   char fosType[LGROUPSIZE][10] = { {0} };
   char fosColour[LGROUPSIZE][20] = { {0} };
   char fosMethod[LGROUPSIZE][20] = { {0} };

   /*  cant get the external reference to 2D char array to work.....
   pickSequences(
       numPS,
       psDepth,
       psVsh,
       psTSF,
       numFOS,
       fosDepth,
       fosType
       ); */

   k = 0;

   for (j = 0; j <= numPS; j++) {
       fprintf(stderr, "PS number = %d \n", j);

       if (j == 0)
       {
           fosDepth[k] = psDepth[j];
           strncpy(fosType[k], "TOP", 3);
           k++;
       }
       else if (j == 1) {
           if (psTSF[j] < psTSF[j - 1] &
               psTSF[j] > psTSF[j + 1])
           {
               fosDepth[k] = psDepth[j];
               strncpy(fosType[k], "HST/TST", 8);
               k++;
           }
           else if (psTSF[j] < psTSF[j - 1] &
               psTSF[j] < psTSF[j + 1])
           {
               fosDepth[k] = psDepth[j];
               strncpy(fosType[k], "HST/HST", 8);
               k++;
           }

       }
       else if (j == numPS )
       {
           fosDepth[k] = psDepth[j];
           strncpy(fosType[k], "BASE", 4);
       }
       else
       {
           fprintf(stderr, "PS number = %d \n", numPS);
           if (psTSF[j] < psTSF[j - 1] &
               psVsh[j - 1] > psVsh[j - 2] * 0.9 &
               psTSF[j - 1] > psTSF[j - 2] &
               psTSF[j] < psTSF[j + 1])
           {

               fosDepth[k] = psDepth[j];
               strncpy(fosType[k], "HST/HST", 8);
               k++;
           }
           else if (psTSF[j] < psTSF[j - 1] &
                    psTSF[j - 1] < psTSF[j - 2] &
                    psVsh[j - 1] < psVsh[j - 2] &
                    psTSF[j] < psTSF[j + 1])
           {

               fosDepth[k] = psDepth[j];
               strncpy(fosType[k], "TST/HST", 8);
               k++;
           }
           else if (psTSF[j] < psTSF[j - 1] &
               psTSF[j - 1] < psTSF[j - 2] &
               psVsh[j - 1] > psVsh[j - 2] &
               psTSF[j] < psTSF[j + 1])
           {

               fosDepth[k] = psDepth[j];
               strncpy(fosType[k], "LST/HST", 8);
               k++;
           }

       }
   }

   numFOS = k;

   // generate thickness, VSH and TSF of each FOS
   calcBedVsh(
       numpts,
       depth,
       vsh,
       numFOS,
       fosDepth,
       fosThick,
       fosVsh,
       fosTSF);

   // generate name of fourth-order sequences
   grpcnt = 1;
   for (k = numFOS - 1; k >= 0; k--, grpcnt++) {
       sprintf(fosName[k], "FOS-%d", grpcnt);
       if (grpcnt % 2 == 0) {
           strncpy(fosColour[k], "yellow_green", 12);
       }
       else {
           strncpy(fosColour[k], "light_sea_green", 15);
       }
   }


   //write out results to .csv file
   FILE* fosfile;
   fosfile = fopen("./data/fos.csv", "w");
   if (fosfile != NULL) {
       fprintf(stderr, "./data/fos.csv has been opened for write \n");
   }

   fprintf(stderr, "number of fourth-order sequences = %d \n", numFOS);

   for (k = 0; k <= numFOS; k++) {
       fprintf(fosfile, "%f, %s, %s, %f, %f, %f, %s \n", fosDepth[k], fosName[k], fosColour[k], fosThick[k], fosVsh[k], fosTSF[k], fosType[k]);
   }

   fclose(fosfile);


}
