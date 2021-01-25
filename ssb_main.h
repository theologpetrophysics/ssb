#ifdef __cplusplus
extern "C" {
#endif

   extern void ssbMain(
        int numpts,
        int halfSmthWindow,
        double *depthLog,
        double *tvdLog,
        double *grLog,
        double *vshLog,
        char* lithgroupMethod,
        double *vshSmth,
        double* vshFirstDeriv,
        double* vshSecondDeriv,
       double* vshSmthFirstDeriv,
       double* vshSmthSecondDeriv,
        double* lithLogValue,
        double lithGroupMinThick); 

    /*
    extern void ssbMain(
        int numpts,
        int halfSmthWindow,
        struct inputLogData inData1,
        double* vshSmth
    ); */

  
#ifdef __cplusplus
}
#endif
