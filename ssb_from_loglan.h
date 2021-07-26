#ifdef __cplusplus
extern "C" {
#endif


    extern void ssbFromLoglan(
        int numpts,
        double depthAveWindow,
        double* depthLog,
        double* tvdLog,
        double* grLog,
        double* vshLog,
        char* lithgroupMethod,
        char* ecsMethod,
        double ecsShaleBreakLimit,
        char* psPickMethod,
        int psIncludeBacksteps,
        double psEsVshRatio,
        double fosNoiseFact,
        char* optElementLog,
        double lithGroupMinThick,
        double* vshSmth,
        double* grSmth,
        double* firstDeriv,
        double* secondDeriv,
        double* lithLogValue);

  
#ifdef __cplusplus
}
#endif
