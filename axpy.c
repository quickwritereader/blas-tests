
#include <bl_test.h>

#undef AXPY

#ifdef COMPLEX
#ifdef DOUBLE
#define AXPY   BLASFUNC(zaxpy)
#else
#define AXPY   BLASFUNC(caxpy)
#endif
#else
#ifdef DOUBLE
#define AXPY   BLASFUNC(daxpy)
#else
#define AXPY   BLASFUNC(saxpy)
#endif
#endif

int main(int argc, char *argv[]) {

    FLOAT *x, *y, *nx, *ny;
    ;
    FLOAT alpha[2] = {2.0, 2.0};
    blasint m, i;
    blasint inc_x = 1, inc_y = 1;
    char *p;

    int zcase = 0;

#ifdef COMPLEX
    //real zero, image zero, both zero, both random
    const char *case_names[] = {"r_zero", "i_zero", "both zero", "both random"};
    int zero_case[] = {1, 2, 3, 0};
#else
    //zero, random
    const char *case_names[] = {"zero", "random"};
    int zero_case[] = {1, 0};
#endif        

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

    LOG( "From : %3d  To : %3d Step = %3d Inc_x = %d Inc_y = %d \n", from, to, step, inc_x, inc_y);

    if ((x = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        LOG( "Out of Memory!!\n");
        exit(1);
    }

    if ((y = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL) {
        LOG( "Out of Memory!!\n");
        exit(1);
    }

    if ((nx = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        LOG( "Out of Memory!!\n");
        exit(1);
    }

    if ((ny = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL) {
        LOG( "Out of Memory!!\n");
        exit(1);
    }

#ifdef linux
    srandom(getpid());
#endif


    for (m = from; m <= to; m += step) {

        for (zcase = 0; zcase<sizeof (zero_case) / sizeof (int); zcase++) {
            int casex = zero_case[zcase];
            LOG( " %6d  %s alpha: \n", (int) m, case_names[zcase]);
            alpha[0] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            alpha[1] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;


            if ((casex & 1) == 1) {
                alpha[0] = 0;
            }
            if ((casex & 2) == 2) {
                alpha[1] = 0;
            }


            for (i = 0; i < m * COMPSIZE * abs(inc_x); i++) {
                x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
                nx[i] = x[i];
            }

            for (i = 0; i < m * COMPSIZE * abs(inc_y); i++) {
                y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
                ny[i] = y[i];
            }


            AXPY(&m, alpha, x, &inc_x, y, &inc_y);
            ref_axpy(m, alpha, nx, inc_x, ny, inc_y);
            compare_vals(m, y, inc_y, ny, inc_y, STRINGIZE(AXPY));

        }
        LOG( "------------\n");

    }
    free(x);
    free(nx);
    free(y);
    free(ny);

    return 0;
}
