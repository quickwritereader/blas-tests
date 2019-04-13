#include <bl_test.h>

#undef GEMM

#ifndef COMPLEX

#ifdef DOUBLE
#define GEMM   BLASFUNC(dgemm)
#else
#define GEMM   BLASFUNC(sgemm)
#endif

#else

#ifdef DOUBLE
#define GEMM   BLASFUNC(zgemm)
#else
#define GEMM   BLASFUNC(cgemm)
#endif

#endif

int main(int argc, char *argv[]) {

    FLOAT *a, *b, *c, *cref, *cbl;
    FLOAT alpha[] = {1.0, 0.0};
    FLOAT beta [] = {1.0, 0.0};
    char transa = 'N';
    char transb = 'N';
    int m, n, k, i, j, lda, ldb, ldc;
    int has_param_m = 0;
    int has_param_n = 0;
    int has_param_k = 0;
    int init_to_one = 0;
    int init_static = 0;
    int lc_add = 0;
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

    if ((p = getenv("BLAS_TRANS"))) {
        transa = *p;
        transb = *p;
    }
    if ((p = getenv("BLAS_TRANSA"))) {
        transa = *p;
    }
    if ((p = getenv("BLAS_TRANSB"))) {
        transb = *p;
    }

    if ((p = getenv("INIT_ONE"))) {
        init_to_one = 1;
    }

    if ((p = getenv("INIT_STATIC"))) {
        init_static = 1;
    }
    if ((p = getenv("BLCABC_ADD"))) {
        lc_add = 1;
    }

    fprintf(stderr, "From : %3d  To : %3d Step=%d : Transa=%c : Transb=%c\n", from, to, step, transa, transb);


    if ((p = getenv("PARAM_M"))) {
        m = atoi(p);
        has_param_m = 1;
    } else {
        m = to;
    }
    if ((p = getenv("PARAM_N"))) {
        n = atoi(p);
        has_param_n = 1;
    } else {
        n = to;
    }
    if ((p = getenv("PARAM_K"))) {
        k = atoi(p);
        has_param_k = 1;
    } else {
        k = to;
    }

    if ((a = (FLOAT *) malloc(sizeof (FLOAT) * (m + lc_add) * (k + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((b = (FLOAT *) malloc(sizeof (FLOAT) * (k + lc_add) * (n + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((c = (FLOAT *) malloc(sizeof (FLOAT) *(m + lc_add) * (n + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((cref = (FLOAT *) malloc(sizeof (FLOAT) * (m + lc_add) * (n + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((cbl = (FLOAT *) malloc(sizeof (FLOAT) * (m + lc_add) * (n + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
#ifdef linux
    srandom(getpid());
#endif
    if (!init_static) {
        for (i = 0; i < (m + lc_add) * (k + lc_add) * COMPSIZE; i++) {
            if (init_to_one) {
                a[i] = 1;
            } else {
                a[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            }
        }
        for (i = 0; i < (k + lc_add) * (n + lc_add) * COMPSIZE; i++) {
            if (init_to_one) {
                b[i] = 1;
            } else {
                b[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            }
        }
        for (i = 0; i < (m + lc_add) * (n + lc_add) * COMPSIZE; i++) {
            if (init_to_one) {
                c[i] = 1;
            } else {
                c[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            }
            cref[i] = c[i];
            cbl[i] = c[i];
        }

    }
    for (i = from; i <= to; i += step) {




        if (!has_param_m) {
            m = i;
        }
        if (!has_param_n) {
            n = i;
        }
        if (!has_param_k) {
            k = i;
        }

        if (transa == 'N') {
            lda = m + lc_add;
        } else {
            lda = k + lc_add;
        }
        if (transb == 'N') {
            ldb = k + lc_add;
        } else {
            ldb = n + lc_add;
        }
        ldc = m + lc_add;

        if (init_static) {
            FLOAT st = 0.0;
            for (i = 0; i < (m + lc_add) * (k + lc_add) * COMPSIZE; i++) {

                if (transa == 'N') {
                    if ((i % lda) == 0) {
                        st = 1.0;
                    }
                    a[i] = st;
                    st += 1.0;
                } else {
                    if ((i % lda) == 0) {
                        st = st + 1.0;
                    }
                    b[i] = st;
                }

            }
            st = 0.0;
            for (i = 0; i < (k + lc_add) * (n + lc_add) * COMPSIZE; i++) {
                if (transb == 'N') {
                    if ((i % ldb) == 0) {
                        st = st + 1.0;
                    }
                    b[i] = st;

                } else {
                    if ((i % ldb) == 0) {
                        st = 1.0;
                    }
                    b[i] = st;
                    st += 1.0;

                }
            }
            for (i = 0; i < (m + lc_add) * (n + lc_add) * COMPSIZE; i++) {

                c[i] = 0;
                cref[i] = c[i];
                cbl[i] = c[i];
            }
        } else {
            for (j = 0; j < (m + lc_add) * (n + lc_add) * COMPSIZE; j++) {
                cref[j] = c[j];
                cbl[j] = c[j];
            }
        }

        fprintf(stderr, " M=%4d, N=%4d, K=%4d lda=%4d ldb=%4d ldc=%4d: \n", (int) m, (int) n, (int) k, lda, ldb, ldc);
        GEMM(&transa, &transb, &m, &n, &k, alpha, a, &lda, b, &ldb, beta, cbl, &ldc);
#ifndef COMPLEX
        ref_gemm(&transa, &transb, &m, &n, &k, alpha, a, &lda, b, &ldb, beta, cref, &ldc);
        compare_vals_aggregated(k, m*n, cbl, 1, cref, 1, STRINGIZE(GEMM));
#endif


    }

    return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
