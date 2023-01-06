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

        bedTSF[j] = bedThick[j] / (1.0 - fmin(bedVsh[j], 0.95));
    }

}




void ssbMain(
    int numpts,
    int halfSmthWindow,
    double* depth,
    double* tvd,
    double* gr,
    double* vsh,
    double* coal,
    int optCoalForce,
    char* lithGroupMethod,
    char* ecsMethod,
    double ecsShaleBreakLimit,
    char* psPickMethod,
    int psIncludeTrendTest,
    double psEsVshRatio,
    double psTsfFactor,
    double fosNoiseFact,
    double fosTsfFactor,
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
    firstDeriv[i] = (derivInLog[i + 1] - derivInLog[i]) / (depth[i + 1] - depth[i]);
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


// First pick sand/shale/coal depths
for (i = 0; i < numpts; i++) {

    //  check for coal
    if (coal[i] == 1.0) {
        lithValue[i] = 2;
    }

    // Force to shale if VSH == 1
    else if (vshSmth[i] == 1.0 & strncmp(optElementLog, "VSH", 3) == 0) {
        lithValue[i] = 1;
    }

    //Force to sand if VSH == 0
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
        else
        {
            // sand lithogroup
            lithValue[i] = 0;
        }
    }
    // use derivative method for lithogroups
    else
    {

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

    else if (lithValue[0] == 2) {
        eValue[j] = 2;
    }

    else {
        eValue[j] = 0;
    }

    j++;

    for (i = 1; i < numpts - 1; i++) {
        if (lithValue[i] != lithValue[i - 1] &
            lithValue[i] == lithValue[i + 1] &
            lithValue[i - 1] == lithValue[i - 2])
        {

            // define current element
            eValue[j] = lithValue[i];
            eDepth[j] = depth[i];

            j++;

            if (debug == 1) {
                fprintf(stderr, "grp depth  = %f \n ", eDepth[j]);
            }

        }
    }

    i = numpts - 1;
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
        else if (eValue[i] == 2) {
            //fprintf(stderr, "lithgroup %d - ", i);
            strncpy(eLith[i], "co", 2);
            strncpy(eColour[i], "black", 6);
        }
        else {
            strncpy(eLith[i], "ss", 2);
            strncpy(eColour[i], "yellow", 6);
        }
    }

    // generate name of elements
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
    double* esType;
    double* esVshRatio;
    double* esTSFRatio;
    allocateMemory1DD(&esDepth, numE, 0);
    allocateMemory1DD(&esVsh, numE, 0);
    allocateMemory1DD(&esThick, numE, 0);
    allocateMemory1DD(&esTSF, numE, 0);
    allocateMemory1DD(&esShaleThick, numE, 0);
    allocateMemory1DD(&esType, numE, 0); // differentiate sand vs coal topped ES
    allocateMemory1DD(&esVshRatio, numE, 0);
    allocateMemory1DD(&esTSFRatio, numE, 0);
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

        else if (i == numE - 1 & eValue[i] == 0) {
            // base surface is base of ES
            esDepth[j] = eDepth[numE];
            esType[j] = 1;
            j++;
        }

        // top of coal is ES
        else if (optCoalForce == 1 & eValue[i] == 2) {
            esDepth[j] = eDepth[i];
            esType[j] = 2;
            j++;
        }

        else if (eValue[i] == 0 &
            eValue[i - 1] == 1 &
            eValue[i + 1] == 1) {
            // top of sand/shale couplet is top of ES unit
            esDepth[j] = eDepth[i];
            esShaleThick[j] = eThick[i + 1];
            esType[j] = 1;
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

    for (k = 1; k <= numES; k++) {
        esVshRatio[k] = esVsh[k - 1] / esVsh[k];
        esTSFRatio[k] = esTSF[k - 1] / esTSF[k];
    }


    //write out results to .csv file
    FILE* esfile;
    esfile = fopen("./data/es.csv", "w");
    if (esfile != NULL) {
        fprintf(stderr, "./data/esfile.csv has been opened for write \n");
    }

    fprintf(stderr, "number of Element Sets = %d \n", numES);

    for (k = 0; k <= numES; k++) {
        fprintf(esfile, "%f, %s, %s, %f, %f, %f, %f, %f, %f \n", esDepth[k], esName[k], esColour[k], esThick[k], esShaleThick[k], esVsh[k], esTSF[k],
            esVshRatio[k], esTSFRatio[k]);
    }

    fclose(esfile);


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  Logic to caluclate ELEMENT COMPLEX SET (ECS) 
    //  We have a choice of Shale break or TSF/VSH methods
    //  Shale break method:
    //    Find when shale is thicker above than below
    //    if shale above current surface is thicker than next shale below then current surface is a shale break
    //  TSF/VSH method:
    //    find where there is an increase in TSF and VSH in overlying Element set
    // 
    //  User has choice of one method or either method


    // top surface is treated as top of ECS
    int numECS;
    int ecsframecnt;
    int estmpcnt;
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
    double* ecsESVshRatio; // Ratio of Vsh of the basal ES to underlying ES for each ECS
    double* ecsESTSFRatio; // Ratio of TSF of the basal ES to underlying ES for each ECS
    double* ecsESNum; // number of ES beds per ECS
    double* ecsType; // indicate is COAL topped ECS
    double* ecsVshRatio;
    double* ecsTSFRatio;

    allocateMemory1DD(&ecsDepth, numES, 0);
    allocateMemory1DD(&ecsVsh, numES, 0);
    allocateMemory1DD(&ecsThick, numES, 0);
    allocateMemory1DD(&ecsSbThick, numES, 0);
    allocateMemory1DD(&ecsSbRatio, numES, 0);
    allocateMemory1DD(&ecsSsThick, numES, 0);
    allocateMemory1DD(&ecsTSF, numES, 0);
    allocateMemory1DD(&ecsESsVsh, numES, 0);
    allocateMemory1DD(&ecsEShVsh, numES, 0);
    allocateMemory1DD(&ecsESVshRatio, numES, 0);
    allocateMemory1DD(&ecsESTSFRatio, numES, 0);
    allocateMemory1DD(&ecsESNum, numES, 0);
    allocateMemory1DD(&ecsType, numES, 0);
    allocateMemory1DD(&ecsVshRatio, numE, 0);
    allocateMemory1DD(&ecsTSFRatio, numE, 0);
    char ecsName[LGROUPSIZE][10] = { {0} };
    char ecsColour[LGROUPSIZE][15] = { {0} };

    k = 0;
    /* ecsthicksum = 0;
     ecsthickmin = 100;
     ecsthickmax = 0;
     ecsvshsum = 0;
   */
    estmpcnt = 0;  //initialise

    for (j = 0; j <= numES; j++) {

        if (j == 0) {
            ecsDepth[k] = esDepth[j];
            k++;
            //fprintf(stderr, "Top at  = %f \n ", ecsDepth[k]);

        }
        else if (j == numES) {
            ecsDepth[k] = esDepth[j];
            //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
            ecsESNum[k-1] = estmpcnt;
            estmpcnt = 0;
            k++;

        }
        // forced to be ECS due to of COAL
        else if (optCoalForce == 1 & esType[j] == 2) {
            ecsDepth[k] = esDepth[j];
            ecsESVshRatio[k] = esVsh[j - 1] / esVsh[j];
            ecsESNum[k - 1] = estmpcnt;
            estmpcnt = 0;
            ecsType[k] = 2.0;
            k++;
        }
        else if ((strncmp(ecsMethod, "SB", 2) == 0 |
            strncmp(ecsMethod, "EITHER", 6) == 0) &
            esShaleThick[j - 1] >= esShaleThick[j] * (1.0 + ecsShaleBreakLimit)) {

            ecsDepth[k] = esDepth[j];
            ecsSbThick[k] = esShaleThick[j - 1];
            ecsSsThick[k] = esThick[j - 1] - esShaleThick[j - 1];
            // calculate shale break thickness change ratio for calibration
            ecsSbRatio[k] = esShaleThick[j - 1] / esShaleThick[j];

            //ecsESVshRatio[k] = esVsh[j - 1] / esVsh[j];
            ecsESVshRatio[k] = esVshRatio[j];
            ecsESTSFRatio[k] = esTSFRatio[j];

            ecsESNum[k-1] = estmpcnt;
            ecsType[k] = 1.0;
            estmpcnt = 0;
            k++;

        }
        else if ((strncmp(ecsMethod, "TSF/VSH", 7) == 0 |
            strncmp(ecsMethod, "EITHER", 6) == 0) &
            esVsh[j - 1] > esVsh[j] &
            esTSF[j - 1] > esTSF[j]) {

            ecsDepth[k] = esDepth[j];
            ecsSbThick[k] = esShaleThick[j - 1];
            ecsSsThick[k] = esThick[j - 1] - esShaleThick[j - 1];
            // calculate shale break thickness change ratio for calibration
            ecsSbRatio[k] = esShaleThick[j - 1] / esShaleThick[j];

            //ecsESVshRatio[k] = esVsh[j - 1] / esVsh[j];
            ecsESVshRatio[k] = esVshRatio[j];
            ecsESTSFRatio[k] = esTSFRatio[j];

            ecsESNum[k-1] = estmpcnt;
            ecsType[k] = 1.0;
            estmpcnt = 0;
            k++;
        }

        estmpcnt++;  // increment the ES counter
    }


    numECS = k - 1;

    // generate name of Element Complex Sets
    grpcnt = 1;
    for (k = numECS - 1; k >= 0; k--, grpcnt++) {
        sprintf(ecsName[k], "ECS-%d", grpcnt);
        //fprintf(stderr, "%s ratio = %f \n ", ecsName[k], ecsESVshRatio[k]);
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

    for (j = 0; j <= numECS; j++) {
        ecsVshRatio[j] = ecsVsh[j - 1] / ecsVsh[j];
        ecsTSFRatio[j] = ecsTSF[j - 1] / ecsTSF[j];
    }

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
        fprintf(ecsfile, "%f, %s, %s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n", ecsDepth[k], ecsName[k], ecsColour[k],
            ecsThick[k], ecsVsh[k], ecsTSF[k], ecsSbThick[k], ecsSbRatio[k], ecsSsThick[k], ecsESsVsh[k], ecsEShVsh[k], ecsESNum[k],
            ecsVshRatio[k], ecsTSFRatio[k], ecsESVshRatio[k], ecsESTSFRatio[k]);
    }

    fclose(ecsfile);


    //-----------------------------------------------------------------------------------
    //  Logic to caluclate PARA-SEQUENCES (PS) 
    //

    int numPS;
    int psframecnt;
    int psecscnt;
    double psthicksum, psthickmin, psthickmax, psvshsum;
    double* psDepth;
    double* psVsh;
    double* psThick;
    double* psTSF;
    double* psECSNum;
    double* psECSVshRatio;
    double* psVshRatio;
    double* psTSFRatio;

    allocateMemory1DD(&psDepth, numECS, 0);
    allocateMemory1DD(&psVsh, numECS, 0);
    allocateMemory1DD(&psThick, numECS, 0);
    allocateMemory1DD(&psTSF, numECS, 0);
    allocateMemory1DD(&psECSNum, numECS, 0);
    allocateMemory1DD(&psECSVshRatio, numECS, 0);
    allocateMemory1DD(&psVshRatio, numECS, 0);
    allocateMemory1DD(&psTSFRatio, numECS, 0);
    char psName[LGROUPSIZE][10] = { {0} };

    char psColour[LGROUPSIZE][20] = { {0} };
    char psMethod[LGROUPSIZE][30] = { {0} };

    k = 0;
    psecscnt = 0;

    for (j = 0; j <= numECS; j++) {
        fprintf(stderr, "%s - ECS TYPE - %f \n", ecsName[j], ecsType[j]);

        if (j == 0) {
            psDepth[k] = ecsDepth[j];
            strncpy(psMethod[k], "Top", 15);
            k++;
            //fprintf(stderr, "Top PS at  = %f \n ", psDepth[k]);

        }
        else if (j == numECS) {
            psDepth[k] = ecsDepth[j];
            strncpy(psMethod[k], "Base", 15);
            psECSNum[k - 1] = psecscnt;
            psecscnt = 0;
            //fprintf(stderr, "Base at  = %f \n ", ecsDepth[k]);
            k++;

        }
        else if (optCoalForce == 1 & ecsType[j] == 2.0) {
            psDepth[k] = ecsDepth[j];
            strncpy(psMethod[k], "Coal", 15);
            psECSNum[k - 1] = psecscnt;
            psecscnt = 0;
            k++;
        }

        // TSF method
        

        else if ((strncmp(psPickMethod, "TSF", 3) == 0 |
            strncmp(psPickMethod, "EITHER", 6) == 0) &
            ecsTSFRatio[j ] > (1.0 + psTsfFactor) &
            ecsVshRatio[j] > (1.0 + psTsfFactor) &
            ecsVsh[j - 1] > 0.30 &
            psecscnt * ecsESNum[j - 1] > 1.0)
             {

          

                psDepth[k] = ecsDepth[j];
                strncpy(psMethod[k], "ECS-TSF/VSH", 30);
                psECSVshRatio[k] = ecsVsh[j - 1] / ecsVsh[j];
                psECSNum[k - 1] = psecscnt;
                psecscnt = 0;
                k++;



        }

        else if ((strncmp(psPickMethod, "TSF", 3) == 0 |
                strncmp(psPickMethod, "EITHER", 6) == 0) &
                ecsTSFRatio[j] >= (1.0) & 
                ecsVshRatio[j] > (1.0 + psTsfFactor) &
                ecsESVshRatio[j] > 1.0 &
                ecsESTSFRatio[j] > 1.0 &
                ecsVsh[j - 1] > 0.30 &
                psecscnt * ecsESNum[j - 1] > 1.0) {

            psDepth[k] = ecsDepth[j];
            strncpy(psMethod[k], "ECS-VSH & ES", 30);
            psECSVshRatio[k] = ecsVsh[j - 1] / ecsVsh[j];
            psECSNum[k - 1] = psecscnt;
            psecscnt = 0;
            k++;
        }

        else if ((strncmp(psPickMethod, "TSF", 3) == 0 |
                strncmp(psPickMethod, "EITHER", 6) == 0) &
                ecsTSFRatio[j] > (1.0 + psTsfFactor) &
                ecsESVshRatio[j] > 1.0 &
                ecsESTSFRatio[j] > 1.0 &
                ecsVsh[j - 1] > 0.30 &
                psecscnt * ecsESNum[j - 1] > 1.0) {

            psDepth[k] = ecsDepth[j];
            strncpy(psMethod[k], "ECS-TSF & ES", 30);
            psECSVshRatio[k] = ecsVsh[j - 1] / ecsVsh[j];
            psECSNum[k - 1] = psecscnt;
            psecscnt = 0;
            k++;
        }

        // shale break method
        else if ((strncmp(psPickMethod, "SB", 2) == 0 |
                strncmp(psPickMethod, "EITHER", 6) == 0) &
                ecsSbThick[j] > ecsSbThick[j + 1] * 1.25 &
                ecsVsh[j - 1] > 0.3 &
                ecsESVshRatio[j] > psEsVshRatio &
                psecscnt * ecsESNum[j - 1] > 1.0)
        {

          

                psDepth[k] = ecsDepth[j];
                strncpy(psMethod[k], "No trend Test - SB/VSH ", 23);
                psECSVshRatio[k] = ecsVsh[j - 1] / ecsVsh[j];
                psECSNum[k - 1] = psecscnt;
                psecscnt = 0;
                k++;


        }

        psecscnt++;
    }

    numPS = k - 1;

    // generate name of Parasequences
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

    for (k = 0; k <= numPS; k++) {
        psVshRatio[k] = psVsh[k - 1] / psVsh[k];
        psTSFRatio[k] = psTSF[k - 1] / psTSF[k];
    }

    //write out results to .csv file
    FILE* psfile;
    psfile = fopen("./data/ps.csv", "w");
    if (psfile != NULL) {
        fprintf(stderr, "./data/ps.csv has been opened for write \n");
    }

    fprintf(stderr, "number of Para-sequences = %d \n", numPS);

    for (k = 0; k <= numPS; k++) {
        fprintf(psfile, "%f, %s, %s, %f, %f, %f, %s, %f, %f, %f, %f \n", psDepth[k], psName[k], psColour[k], psThick[k], psVsh[k], psTSF[k], 
            psMethod[k], psECSNum[k], psVshRatio[k], psTSFRatio[k], psECSVshRatio[k]);
    }

    fclose(psfile);



    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //  Now need to find fourth order surface tops....
    //  

    int numFOS, numFOStmp, currentECSIndex, l;
    double* fosDepth;
    double* fosDepthtmp;
    double* fosVsh;
    double* fosThick;
    double* fosTSF;
    double* fosDel;
    double* sVshRatio;
    double* sTSFRatio;
    double lastMtsPSVsh, vshTmp, lastFOS; // lastFOS 1 = mts, 0 = mrs
    allocateMemory1DD(&fosDepth, numPS, 0);
    allocateMemory1DD(&fosDepthtmp, numPS, 0);
    allocateMemory1DD(&fosVsh, numPS, 0);
    allocateMemory1DD(&fosThick, numPS, 0);
    allocateMemory1DD(&fosTSF, numPS, 0);
    allocateMemory1DD(&fosDel, numPS, 0);
    allocateMemory1DD(&sVshRatio, numPS, 0);
    allocateMemory1DD(&sTSFRatio, numPS, 0);
    char fosName[LGROUPSIZE][10] = { {0} };
    char fosType[LGROUPSIZE][12] = { {0} };
    char fosTypetmp[LGROUPSIZE][10] = { {0} };
    char fosColour[LGROUPSIZE][20] = { {0} };
    char fosMethod[LGROUPSIZE][20] = { {0} };
    char fosSysTract[LGROUPSIZE][20] = { {0} };
    //int fosTsfScore, fosVshScore, fosEcsScore;

    
    k = 0;
    lastMtsPSVsh = 0;

    // first we look for the FOS surfaces using
    //  PS TSF and ECS Vsh
    // need to refer to the ECS vsh ratio at PS boundary previously generated - psECSVshRatio
    for (j = 0; j <= numPS; j++) {


        //find the equivalent ECS index
        for (l = 0; l < numECS; l++) {
            if (ecsDepth[l] == psDepth[j]) {
                currentECSIndex = l;
            }
        }
        fprintf(stderr, "PS number = %s \n", psName[j]);
        fprintf(stderr, "ECS number = %s \n", ecsName[currentECSIndex]);

        //Deal with top and basal PS units
        if (j == 0)
        {
            fosDepthtmp[k] = psDepth[j];
            strncpy(fosTypetmp[k], "TOP", 4);
            k++;
        }
        else if (j == numPS)
        {
            fosDepthtmp[k] = psDepth[j];
            strncpy(fosTypetmp[k], "BASE", 4);
            k++;
        }
        /*else if (psTSF[j] * (1.0 + fosTsfFactor) < psTSF[j - 1] &
                psECSVshRatio[j] > 0.975 &
                psVsh[j - 1] * (1.0 + fosNoiseFact) > psVsh[j])*/

        else if (psTSFRatio[j] > (1.0 + fosTsfFactor) &
                psVshRatio[j] > (1.0 - fosNoiseFact) &
                psECSVshRatio[j] > (1.0 - fosNoiseFact))

        {
            fosDepthtmp[k] = psDepth[j];
            strncpy(fosTypetmp[k], "Mts ", 4);
            vshTmp = lastMtsPSVsh;
            lastMtsPSVsh = psVsh[j];
            lastFOS = 1;
            // tag for deletion if directly below another one
            if (psDepth[j - 1] == fosDepthtmp[k - 1])
            {
                fosDel[k] = 1.0;
                lastMtsPSVsh = vshTmp;
            }

            k++;

            fprintf(stderr, "FOS - PS number = %s \n", psName[j]);
        }
        // find the Mrs surfaces
        // increase in ECStsf and ECSVsh
        // PS Vsh below Mts > PSVsh below Mrs
         else if (ecsTSF[currentECSIndex] * 1.5 < ecsTSF[currentECSIndex-1] &
                  ecsVsh[currentECSIndex] < ecsVsh[currentECSIndex - 1] &
                  psVsh[j] < lastMtsPSVsh &
                  lastFOS == 1)
        {
            
            fosDepthtmp[k] = psDepth[j];
            strncpy(fosTypetmp[k], "Mrs ", 4);
            lastFOS = 0;

            k++;
        }
     }
    numFOStmp = k;


    //delete the repeated ones

    k = 0;

    for (j = 0; j < numFOStmp; j++) {

        if (fosDel[j] != 1.0)
        {
            fosDepth[k] = fosDepthtmp[j];
            strncpy(fosType[k], fosTypetmp[j], 4);
            k++;
        }
    }
    numFOS = k-1;




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

    // generate name of fourth-order surfaces and system tracks
    grpcnt = 1;
    for (k = numFOS - 1; k >= 0; k--, grpcnt++) {
        sprintf(fosName[k], "FOS-%d", grpcnt);

        
        if (strncmp(fosType[k + 1], "BASE", 4) == 0)
        {
            if (strncmp(fosType[k], "Mts", 3) == 0)
            {
                strncpy(fosColour[k], "DARK_GREEN", 15);
                strncpy(fosSysTract[k], "TST", 12);
            } 
            
            else if (strncmp(fosType[k], "Mrs", 3) == 0)
            {
                strncpy(fosColour[k], "BLUE", 15);
                strncpy(fosSysTract[k], "RST", 12);
            }

            else {
                strncpy(fosColour[k], "GRAY", 15);
                strncpy(fosSysTract[k], "?", 12);
            }
        } 
        
        else if (strncmp(fosType[k], "TOP", 3) == 0)
        {
            if (strncmp(fosType[k + 1], "Mts", 3) == 0)
            {
                strncpy(fosColour[k], "BLUE", 15);
                strncpy(fosSysTract[k], "RST", 12);
            }
            else if (strncmp(fosType[k], "Mrs", 3) == 0)
            {
                strncpy(fosColour[k], "DARK_GREEN", 15);
                strncpy(fosSysTract[k], "TST", 12);
            }
        }
           
        else if (strncmp(fosType[k], "Mts", 3) == 0 &
                 strncmp(fosType[k+1], "Mrs", 3) == 0)
        {
            strncpy(fosColour[k], "dark_green", 12);
            strncpy(fosSysTract[k], "TST", 12);
        }

        else if (strncmp(fosType[k], "Mrs", 3) == 0 &
                 strncmp(fosType[k + 1], "Mts", 3) == 0)
        {
            strncpy(fosColour[k], "blue", 15);
            strncpy(fosSysTract[k], "RST", 12);
        }

        else if (strncmp(fosType[k], "Mts", 3) == 0 &
                 strncmp(fosType[k + 1], "Mts", 3) == 0)
        {
            strncpy(fosColour[k], "blue", 15);
            strncpy(fosSysTract[k], "RST", 12);
        }

        else 
         {
            strncpy(fosColour[k], "GRAY", 15);
            strncpy(fosSysTract[k], "??", 12);
         }

    }


    //write out results to .csv file
    FILE* fosfile;
    fosfile = fopen("./data/fos.csv", "w");
    if (fosfile != NULL) {
        fprintf(stderr, "./data/fos.csv has been opened for write \n");
    }

    fprintf(stderr, "number of fourth-order surfaces = %d \n", numFOS);

    for (k = 0; k <= numFOS; k++) {
        fprintf(fosfile, "%f, %s, %s, %f, %f, %f, %s, %s \n", fosDepth[k], fosName[k], fosColour[k], fosThick[k], 
            fosVsh[k], fosTSF[k], fosType[k], fosSysTract[k]);
    }

    fclose(fosfile);


    // generate the sequences
    int numS, currentPSIndex;
    double* sDepth;
    double* sVsh;
    double* sThick;
    double* sTSF;
    double* sPsVshRatio;
    allocateMemory1DD(&sDepth, numFOS, 0);
    allocateMemory1DD(&sVsh, numFOS, 0);
    allocateMemory1DD(&sThick, numFOS, 0);
    allocateMemory1DD(&sTSF, numFOS, 0);
    allocateMemory1DD(&sPsVshRatio, numFOS, 0);
    char sName[LGROUPSIZE][10] = { {0} };
    char sType[LGROUPSIZE][10] = { {0} };
    char sColour[LGROUPSIZE][20] = { {0} };

    k = 0;

    for (j = 0; j <= numFOS; j++) {

        if (j == 0) {
            sDepth[k] = fosDepth[j];
            strncpy(sType[k], "TOP", 3);
            k++;
        }
        else if (j == numFOS) {
            sDepth[k] = fosDepth[j];
            strncpy(sType[k], "BASE", 4);
            k++;
        }
        else if (strncmp(fosType[j], "Mrs", 3) == 0 |
            (strncmp(fosType[j], "Mts", 3) == 0 &
            strncmp(fosType[j + 1], "Mts", 3) == 0) ) 
        {
            sDepth[k] = fosDepth[j];
            strncpy(sType[k], "Sb", 9);
            k++;
        }
    }

    numS = k - 1;

    // generate thickness, VSH and TSF of each Sequence
    calcBedVsh(
        numpts,
        depth,
        vsh,
        numS,
        sDepth,
        sThick,
        sVsh,
        sTSF);

    for (k = 0; k <= numS; k++) {
        sVshRatio[k] = sVsh[k - 1] / sVsh[k];
        sTSFRatio[k] = sTSF[k - 1] / sTSF[k];

        //find the equivalent PS index
        for (l = 0; l < numPS; l++) {
            if (psDepth[l] == sDepth[j]) {
                currentPSIndex = l;
            }
        }

        sPsVshRatio[k] = psVsh[currentPSIndex - 1] / psVsh[currentPSIndex];

    }

    // generate name of fourth-order sequences
    grpcnt = 1;
    for (k = numS - 1; k >= 0; k--, grpcnt++) {
        sprintf(sName[k], "S-%d", grpcnt);
        if (grpcnt % 2 == 0) {
            strncpy(sColour[k], "dark_sea_green", 18);
        }
        else {
            strncpy(sColour[k], "light_goldenrod", 18);
        }
    }

    //write out results to .csv file
    FILE* sfile;
    sfile = fopen("./data/s.csv", "w");
    if (sfile != NULL) {
        fprintf(stderr, "./data/s.csv has been opened for write \n");
    }

    fprintf(stderr, "number of sequences = %d \n", numS);

    for (k = 0; k <= numS; k++) {
        fprintf(sfile, "%f, %s, %s, %f, %f, %f, %s, %f \n", sDepth[k], sName[k], sColour[k], sThick[k], sVsh[k], sTSF[k], sType[k], sPsVshRatio[k]);
    }

    fclose(sfile);


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //  Now need to find third order sequences tops....
      //  

    int numTOS;
    double* tosDepth;
    double* tosVsh;
    double* tosThick;
    double* tosTSF;
    double* tosDel;
    //double lastMtsFOSVsh, lastTOS; // lastTOS 1 = mts, 0 = mrs
    double tosTsfFactor, tosNoiseFact;
    allocateMemory1DD(&tosDepth, numS, 0);
    allocateMemory1DD(&tosVsh, numS, 0);
    allocateMemory1DD(&tosThick, numS, 0);
    allocateMemory1DD(&tosTSF, numS, 0);
    allocateMemory1DD(&tosDel, numS, 0);
    char tosName[LGROUPSIZE][10] = { {0} };
    char tosType[LGROUPSIZE][12] = { {0} };
    char tosColour[LGROUPSIZE][20] = { {0} };
    char tosMethod[LGROUPSIZE][20] = { {0} };
    char tosSysTract[LGROUPSIZE][20] = { {0} };

    k = 0;

    tosTsfFactor = 0.25;
    tosNoiseFact = 0.0;

    // first we look for the TOS surfaces using
    //  Seq TSF and ECS Vsh
    //  may need to refer to the PS vsh ratio at S boundary previously generated - sPSVshRatio

    for (j = 0; j <= numS; j++) {


        //find the equivalent PS index
        for (l = 0; l < numPS; l++) {
            if (psDepth[l] == sDepth[j]) {
                currentPSIndex = l;
            }
        }
        
        
        fprintf(stderr, "j number = %d \n", j);
        fprintf(stderr, "PS number = %s \n", psName[currentPSIndex]);

        //Deal with top and basal PS units
        if (j == 0)
        {
            tosDepth[k] = sDepth[j];
            strncpy(tosType[k], "TOP", 4);
            k++;
        }
        else if (j == numS)
        {
            tosDepth[k] = sDepth[j];
            strncpy(tosType[k], "BASE", 4);
            fprintf(stderr, "K @ base = %d\n", k);
            k++;
        }

        else if (sTSFRatio[j] > (1.0 + tosTsfFactor) &
            (sVshRatio[j] > (1.0 - tosNoiseFact) |
            sPsVshRatio[j] > (1.0 - tosNoiseFact)))

        {
            if (sVsh[j] > sVsh[j - 1])
            {
                tosDepth[k] = sDepth[j];
                strncpy(tosType[k], "MTS/MRS ", 8);
                k++;

            }
            else {

                tosDepth[k] = sDepth[j];
                strncpy(tosType[k], "MRS ", 4);
                k++;
            }


        }
    }

    numTOS = k-1;

    // generate thickness, VSH and TSF of each TOS
    calcBedVsh(
        numpts,
        depth,
        vsh,
        numTOS,
        tosDepth,
        tosThick,
        tosVsh,
        tosTSF);

    // generate name of third-order sequences
    grpcnt = 1;
    for (k = numTOS - 1; k >= 0; k--, grpcnt++) {
        sprintf(tosName[k], "TOS-%d", grpcnt);


        if (strncmp(tosType[k + 1], "BASE", 4) == 0)
        {
            strncpy(tosColour[k], "royal_blue", 15);
            strncpy(tosSysTract[k], "RSS", 12);
        }

        else if (strncmp(tosType[k], "TOP", 3) == 0)
        {
            if (strncmp(tosType[k + 1], "MTS", 3) == 0)
            {
                strncpy(tosColour[k], "royal_blue", 15);
                strncpy(tosSysTract[k], "RSS", 12);
            }
            else if (strncmp(tosType[k+1], "MRS", 3) == 0)
            {
                strncpy(tosColour[k], "sea_green", 15);
                strncpy(tosSysTract[k], "TSS", 12);
            }
        }

        else if (strncmp(tosType[k], "MTS", 3) == 0 &
            strncmp(tosType[k + 1], "MRS", 3) == 0)
        {
            strncpy(tosColour[k], "sea_green", 12);
            strncpy(tosSysTract[k], "RSS", 12);
        }

        else if (strncmp(tosType[k], "MRS", 3) == 0 &
            strncmp(tosType[k + 1], "MTS", 3) == 0)
        {
            strncpy(tosColour[k], "royal_blue", 15);
            strncpy(tosSysTract[k], "RSS", 12);
        }

        else if (strncmp(tosType[k], "MTS", 3) == 0 &
            strncmp(tosType[k + 1], "MTS", 3) == 0)
        {
            strncpy(tosColour[k], "royal_blue", 15);
            strncpy(tosSysTract[k], "RSS", 12);
        }

        else
        {
            strncpy(tosColour[k], "GRAY", 15);
            strncpy(tosSysTract[k], "??", 12);
        }

    }


    //write out results to .csv file
    FILE* tosfile;
    fprintf(stderr, "./data/tos.csv trying to open \n");
    tosfile = fopen("./data/tos.csv", "w");
    if (tosfile != NULL) {
        fprintf(stderr, "./data/tos.csv has been opened for write \n");
    }

    fprintf(stderr, "number of third-order surfaces = %d \n", numTOS);

    for (k = 0; k <= numTOS; k++) {
        fprintf(tosfile, "%f, %s, %s, %f, %f, %f, %s, %s \n", tosDepth[k], tosName[k], tosColour[k], tosThick[k],
            tosVsh[k], tosTSF[k], tosType[k], tosSysTract[k]);
    }

    fclose(tosfile);
    fprintf(stderr, "./data/tos.csv closed \n");


    // generate the sequence sets
    int numSS;
    double* ssDepth;
    double* ssVsh;
    double* ssThick;
    double* ssTSF;
    allocateMemory1DD(&ssDepth, numTOS, 0);
    allocateMemory1DD(&ssVsh, numTOS, 0);
    allocateMemory1DD(&ssThick, numTOS, 0);
    allocateMemory1DD(&ssTSF, numTOS, 0);
    char ssName[LGROUPSIZE][10] = { {0} };
    char ssType[LGROUPSIZE][10] = { {0} };
    char ssColour[LGROUPSIZE][20] = { {0} };

    k = 0;

    for (j = 0; j <= numTOS; j++) {

        if (j == 0) {
            ssDepth[k] = tosDepth[j];
            strncpy(ssType[k], "TOP", 3);
            k++;
        }
        else if (j == numTOS) {
            ssDepth[k] = tosDepth[j];
            strncpy(ssType[k], "BASE", 4);
            k++;
        }
        else if (strncmp(tosType[j], "MRS", 3) == 0)
        {
            ssDepth[k] = tosDepth[j];
            strncpy(ssType[k], "Sb", 9);
            k++;
        }
        else if (strncmp(tosType[j], "MTS", 3) == 0 &
                strncmp(tosType[j + 1], "MTS", 3) == 0)
        {
            ssDepth[k] = tosDepth[j];
            strncpy(ssType[k], "Sb", 9);
            k++;
        }
        else if (strncmp(tosType[j], "MTS", 3) == 0 &
            strncmp(tosType[j + 1], "BASE", 4) == 0)
            {
                ssDepth[k] = tosDepth[j];
                strncpy(ssType[k], "Sb", 9);
                k++;
        }
    }

    numSS = k - 1;

    // generate thickness, VSH and TSF of each Sequence
    calcBedVsh(
        numpts,
        depth,
        vsh,
        numSS,
        ssDepth,
        ssThick,
        ssVsh,
        ssTSF);

    // generate name of third-order sequences
    grpcnt = 1;

    if (numSS > 1)
        {
            for (k = numSS - 1; k >= 0; k--, grpcnt++) {
                sprintf(ssName[k], "SS-%d", grpcnt);
                if (grpcnt % 2 == 0) {
                    strncpy(ssColour[k], "light_steel_blue", 18);
                }
                else {
                    strncpy(ssColour[k], "royal_blue1", 18);
                }
            }

            //write out results to .csv file
            FILE* ssfile;
            ssfile = fopen("./data/ss.csv", "w");
            if (ssfile != NULL) {
                fprintf(stderr, "./data/ss.csv has been opened for write \n");
            }

            fprintf(stderr, "number of sequence sets = %d \n", numSS);

            for (k = 0; k <= numSS; k++) {
                fprintf(ssfile, "%f, %s, %s, %f, %f, %f, %s \n", ssDepth[k], ssName[k], ssColour[k], ssThick[k], ssVsh[k], ssTSF[k], ssType[k]);
            }

            fclose(ssfile);
        }


}