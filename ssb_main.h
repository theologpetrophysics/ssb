#ifdef __cplusplus
extern "C" {
#endif

   extern void ssbMain(
        int numpts,
        int halfSmthWindow,
        double* depthLog,
        double* tvdLog,
        double* grLog,
        double* vshLog,
        char* lithgroupMethod,
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
