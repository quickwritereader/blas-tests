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
    int m, n, k, i, j;
    int has_param_m = 0;
    int has_param_n = 0;
    int has_param_k = 0;
    int init_to_one = 0;
    int init_static = 0;
    int ldc=1;
    int ldb=1;
    int lda=1;
    char *p;

    int from = 1;
    int to = 1;
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
    if ((p = getenv("BLAS_LDA"))) {
        lda = atoi(p);
    }
    if ((p = getenv("BLAS_LDB"))) {
        ldb =  atoi(p);
    }
     if ((p = getenv("BLAS_LDC"))) {
        ldc =  atoi(p);
    }
 

   if((p=getenv("ZERO_BETA"))){
             beta[0]=0;
	     beta[1]=0;
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
    unsigned char nota, notb;
    int nrowa,ncola, nrowb, ncolb;
    nota = (transa == 'N');
    notb = (transb == 'N');
    if (nota) {
        nrowa = m;
        ncola = k;
    } else {
        nrowa = k;
        ncola = m;
    }
    if (notb) {
        nrowb = k;
        ncolb = n;
    } else {
        nrowb = n;
        ncolb = k;
    }

    lda = MAX(lda, MAX(1,nrowa));
    ldb = MAX(ldb, MAX(1,nrowb));
    ldc = MAX(ldc, MAX(1,m));
    fprintf(stdout,"%d %d %d with leading : %d %d %d\n",m,n,k,lda,ldb,ldc);

    if ((a = (FLOAT *) malloc(sizeof (FLOAT) * lda * ncola * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((b = (FLOAT *) malloc(sizeof (FLOAT) * ldb * ncolb * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((c = (FLOAT *) malloc(sizeof (FLOAT) * ldc * n * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((cref = (FLOAT *) malloc(sizeof (FLOAT) * ldc * n * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((cbl = (FLOAT *) malloc(sizeof (FLOAT) * ldc * n * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }    
#ifdef linux
    srandom(getpid());
#endif
    if (!init_static) {
        for (i = 0; i < lda * ncola * COMPSIZE; i++) {
            if (init_to_one) {
                a[i] = 1;
            } else {
                a[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            }
        }
        for (i = 0; i < ldb * ncolb * COMPSIZE; i++) {
            if (init_to_one) {
                b[i] = 1;
            } else {
                b[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            }
        }
        for (i = 0; i < ldc * n * COMPSIZE; i++) {
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

 

        if (init_static) {
            FLOAT st = 0.0;
            for (i = 0; i < lda * ncola * COMPSIZE; i++) {

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
                    a[i] = st;
                }

            }
            st = 0.0;
            for (i = 0; i < ldb * ncolb * COMPSIZE; i++) {
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
            for (i = 0; i < ldc * n * COMPSIZE; i++) {
                c[i] = 0;
                cref[i] = c[i];
                cbl[i] = c[i];
            }
        } else {
            for (j = 0; j < ldc * n * COMPSIZE; j++) {
                cref[j] = c[j];
                cbl[j] = c[j];
            }
        }

        fprintf(stderr, " M=%4d, N=%4d, K=%4d lda=%4d ldb=%4d ldc=%4d: \n", (int) m, (int) n, (int) k, lda, ldb, ldc);
        GEMM(&transa, &transb, &m, &n, &k, alpha, a, &lda, b, &ldb, beta, cbl, &ldc);
#ifndef COMPLEX
        ref_gemm(&transa, &transb, &m, &n, &k, alpha, a, &lda, b, &ldb, beta, cref, &ldc); 
#else
         ref_zgemm(&transa, &transb, &m, &n, &k,(COMPLEX_FLOAT*) alpha, (COMPLEX_FLOAT*) a, &lda, (COMPLEX_FLOAT*) b, &ldb, (COMPLEX_FLOAT*)beta, cref, &ldc);
                    
#endif
         compare_vals_aggregated(k, m*n, cbl, 1, cref, 1, STRINGIZE(GEMM));  

    }
    free(c);
    free(a);
    free(b);
    free(cbl);
    free(cref);
    return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
