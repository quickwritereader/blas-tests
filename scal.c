
#include <bl_test.h>

#undef SCAL

#ifdef COMPLEX
#ifdef DOUBLE
#define SCAL   BLASFUNC(zscal)
#else
#define SCAL   BLASFUNC(cscal)
#endif
#else
#ifdef DOUBLE
#define SCAL   BLASFUNC(dscal)
#else
#define SCAL   BLASFUNC(sscal)
#endif
#endif

int main(int argc, char *argv[]) {

    FLOAT *x, *nx;
    ;
    FLOAT alpha[2] = {2.0, 2.0};
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

    blasint m, i;
    blasint inc_x = 1;
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


    fprintf(stderr, "From : %3d  To : %3d Step = %3d Inc_x = %d  \n", from, to, step, inc_x);

    if ((x = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }


    if ((nx = (FLOAT *) malloc(sizeof (FLOAT) * to * abs(inc_x) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }



#ifdef linux
    srandom(getpid());
#endif


    for (m = from; m <= to; m += step) {

        for (zcase = 0; zcase<sizeof (zero_case) / sizeof (int); zcase++) {
            int casex = zero_case[zcase];
            fprintf(stderr, " %6d  %s alpha: \n", (int) m, case_names[zcase]);
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




            SCAL(&m, alpha, x, &inc_x);

            ref_scal(m, alpha, nx, inc_x);
            compare_vals(m, x, inc_x, nx, inc_x, "SCAL: ");


        }
        fprintf(stderr, "------------\n");

    }

    free(x);
    free(nx);
    return 0;
}




