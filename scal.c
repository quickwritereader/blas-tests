
#include <bl_test.h>

#undef SCAL


#ifdef DOUBLE
#define SCAL   BLASFUNC(dscal)
#else
#define SCAL   BLASFUNC(sscal)
#endif

int main(int argc, char *argv[]) {

    FLOAT *x,   *nx ;
    ;
    FLOAT alpha[2] = {2.0, 2.0};
    blasint m, i;
    blasint inc_x = 1 ;
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
 

    fprintf(stderr, "From : %3d  To : %3d Step = %3d Inc_x = %d  \n", from, to, step, inc_x );

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

        fprintf(stderr, " %6d : ", (int) m);


        for (i = 0; i < m * COMPSIZE * abs(inc_x); i++) {
            x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            nx[i] = x[i];
        }
 

        alpha[0] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        alpha[1] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;

        SCAL (&m, alpha, x, &inc_x);
        ref_scal(m, alpha[0], nx, inc_x);
        compare_vals(m, x,inc_x, nx,inc_x, "SCAL: ");
        fprintf(stderr, "------------\n");

    }
 

return 0;
}
 

    	

    	 