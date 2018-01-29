
#include <bl_test.h>

#undef COPY

#ifdef COMPLEX
#ifdef DOUBLE
#define COPY   BLASFUNC(zcopy)
#else
#define COPY   BLASFUNC(ccopy)
#endif
#else
#ifdef DOUBLE
#define COPY   BLASFUNC(dcopy)
#else
#define COPY   BLASFUNC(scopy)
#endif
#endif

int main(int argc, char *argv[]) {

    FLOAT *x, *y;

    blasint m, i;
    blasint inc_x = 1;
    blasint inc_y = 1;
    char *p;

    int from = 1;
    int to = 200;
    int step = 1;


    argc--;
    argv++;

    if (argc > 0) {
        from = atol(*argv);
        argc--;
        argv++;
    }
    if (argc > 0) {
        to = MAX(atol(*argv), from);
        argc--;
        argv++;
    }
    if (argc > 0) {
        step = atol(*argv);
        argc--;
        argv++;
    }


    if ((p = getenv("OPENBLAS_INCX"))) inc_x = atoi(p);
    if ((p = getenv("OPENBLAS_INCY"))) inc_y = atoi(p);

    LOG( "From : %3d  To : %3d Step = %3d Inc_x = %d  Inc_y = %d \n", from, to, step, inc_x, inc_y);

    if ((x = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        LOG( "Out of Memory!!\n");
        exit(1);
    }


    if ((y = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL) {
        LOG( "Out of Memory!!\n");
        exit(1);
    }



#ifdef linux
    srandom(getpid());
#endif


    for (m = from; m <= to; m += step) {

        LOG( " %6d : ", (int) m);


        for (i = 0; i < m * COMPSIZE * abs(inc_x); i++) {
            x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        }

        for (i = 0; i < m * COMPSIZE * abs(inc_y); i++) {
            y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        }



        COPY(&m, x, &inc_x, y, &inc_y);

        compare_vals(m, x, inc_x, y, inc_y, STRINGIZE(COPY));

        LOG( "------------\n");

    }

    free(x);
    free(y);
    return 0;
}




