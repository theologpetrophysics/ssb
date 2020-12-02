#ifdef __cplusplus
extern "C" {
#endif

    extern void ssbFromLoglan(
        int numpts,
        int smthFactor,
        double *depthLog,
        double *tvdLog,
        double *grLog,
        double *vshLog,
        double *vshSmth
    );

    extern int smoothLogData(
        int numFrames,
        char* wtShape,
        double* logData,
        double* sValue);

  
#ifdef __cplusplus
}
#endif
