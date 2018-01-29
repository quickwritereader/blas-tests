#include "bl_test.h"
#include <complex.h>
#undef DOT
#ifdef COMPLEX
#ifdef CONJ
#ifdef DOUBLE
#define DOT   BLASFUNC(zdotc) 
#else
#define DOT   BLASFUNC(cdotc) 
#endif
#else
#ifdef DOUBLE
#define DOT   BLASFUNC(zdotu) 
#else
#define DOT   BLASFUNC(cdotu) 
#endif
#endif


#else
#ifdef DOUBLE
#define DOT   BLASFUNC(ddot) 
#else
#define DOT   BLASFUNC(sdot) 
#endif
#endif

int main(int argc, char *argv[]) {

    FLOAT *x, *y;

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

    LOG( "From : %3d  To : %3d Step = %3d Inc_x = %d Inc_y = %d \n", from, to, step, inc_x, inc_y);

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


        LOG( " %6d : \n", (int) m);


        for (i = 0; i < m * COMPSIZE * abs(inc_x); i++) {
            x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        }

        for (i = 0; i < m * COMPSIZE * abs(inc_y); i++) {
            y[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
        }

#ifdef COMPLEX
        FLOAT _Complex result=DOT(&m, x, &inc_x, y, &inc_y);
        COMPLEX_FLOAT result2 = ref_zcdot(m, x, inc_x, y, inc_y); 
#else        
        FLOAT result = DOT(&m, x, &inc_x, y, &inc_y);
        FLOAT result2 = ref_dot(m, x, inc_x, y, inc_y);
#endif
       compare_aggregate(m, (FLOAT*)&result,  (FLOAT*)&result2, STRINGIZE(DOT));
        LOG( "------------\n");


    }
    free(x);
    free(y);
    return 0;
}


