#include "bl_test.h"

#undef DOT

#ifdef DOUBLE
#define DOT   BLASFUNC(ddot) 
#else
#define DOT   BLASFUNC(sdot) 
#endif

int main(int argc, char *argv[]) {

    FLOAT *x, *y;
    FLOAT result;
    blasint m, i;
    blasint inc_x = 1, inc_y = 1;


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

    fprintf(stderr, "From : %3d  To : %3d Step = %3d Inc_x = %d Inc_y = %d \n", from, to, step, inc_x, inc_y);

    if ((x = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }

    if ((y = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }

#ifdef linux
    srandom(getpid());
#endif



    for (m = from; m <= to; m += step) {


        fprintf(stderr, " %6d : \n", (int) m);


        for (i = 0; i < m * COMPSIZE * abs(inc_x); i++) {
            x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        }

        for (i = 0; i < m * COMPSIZE * abs(inc_y); i++) {
            y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        }


        result = DOT(&m, x, &inc_x, y, &inc_y);
        FLOAT result2 = ref_dot(m, x, inc_x, y, inc_y);

        compare_vals(1, &result,1, &result2,1, "DOT ERR: ");
        fprintf(stderr, "------------\n");


    }

    return 0;
}


