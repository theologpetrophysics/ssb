#ifdef __cplusplus
extern "C" {
#endif


    extern void ssbFromLoglan(
        int numpts,
        double depthAveWindow,
        double *depthLog,
        double *tvdLog,
        double *grLog,
        double *vshLog,
        char* lithgroupMethod,
        double lithGroupMinThick,
        double *vshSmth,
        double* vshFirstDeriv,
        double* vshSecondDeriv,
        double* vshSmthFirstDeriv,
        double* vshSmthSecondDeriv,
        double *lithLogValue
    );

  
#ifdef __cplusplus
}
#endif
