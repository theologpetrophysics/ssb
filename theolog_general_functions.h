#ifdef __cplusplus
extern "C" {
#endif

    extern double smoothLogData(
        int numFrames,
        char* wtShape,
        double* logData);

    extern double smoothLogDataSimple(
        int numFrames,
        double* logData);

    extern int allocateMemory1DD(
        double** ptr,
        int n,
        int feedback);

    extern int allocateMemory1DI(
        int** ptr,
        int n,
        int feedback);

  
#ifdef __cplusplus
}
#endif
