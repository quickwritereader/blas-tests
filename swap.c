
#include <bl_test.h>

#undef SWAP

#ifdef COMPLEX
#ifdef DOUBLE
#define SWAP   BLASFUNC(zswap)
#else
#define SWAP   BLASFUNC(cswap)
#endif
#else
#ifdef DOUBLE
#define SWAP   BLASFUNC(dswap)
#else
#define SWAP   BLASFUNC(sswap)
#endif
#endif

int main(int argc, char *argv[]) {

    FLOAT *x, *y, *nx, *ny;

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

    fprintf(stderr, "From : %3d  To : %3d Step = %3d Inc_x = %d  Inc_y = %d \n", from, to, step, inc_x, inc_y);

    if ((x = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }


    if ((y = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((nx = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }


    if ((ny = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_y) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }


#ifdef linux
    srandom(getpid());
#endif


    for (m = from; m <= to; m += step) {

        fprintf(stderr, " %6d : ", (int) m);


        for (i = 0; i < m * COMPSIZE * abs(inc_x); i++) {
            x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            nx[i] = x[i];
        }

        for (i = 0; i < m * COMPSIZE * abs(inc_y); i++) {
            y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            ny[i] = y[i];
        }



        SWAP(&m, x, &inc_x, y, &inc_y);

        compare_vals(m, x, inc_x, ny, inc_y, "SWAP X: ");
        compare_vals(m, nx, inc_x, y, inc_y, "SWAP Y: ");
        fprintf(stderr, "------------\n");

    }
    free(x);
    free(nx);
    free(y);
    free(ny);

    return 0;
}




