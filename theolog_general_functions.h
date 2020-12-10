#ifdef __cplusplus
extern "C" {
#endif

  /*  extern double smoothLogData(
        int numFrames,
        char* wtShape,
        double* logData); */

    extern double smoothLogDataSimple(
        int numFrames,
        double* logData);

    extern int allocateMemory(
        double** ptr,
        int n);

  
#ifdef __cplusplus
}
#endif
