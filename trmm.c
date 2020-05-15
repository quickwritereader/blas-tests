#include <bl_test.h>

#undef TRMM

#ifndef COMPLEX

#ifdef DOUBLE
#define TRMM BLASFUNC(dtrmm)
#else
#define TRMM BLASFUNC(strmm)
#endif

#else

#ifdef DOUBLE
#define TRMM BLASFUNC(ztrmm)
#else
#define TRMM BLASFUNC(ctrmm)
#endif

#endif

int main(int argc, char* argv[]) {
  FLOAT *a, *b, *bi,*bref;
  FLOAT alpha[] = {1.0, 1.0}; 
  int has_param_n = 0;
  int has_param_m = 0;
  int init_to_one = 0;
  int init_static = 0;
  int lc_add = 0;
  char* p;

  char side = 'L';
  char uplo = 'U';
  char trans = 'N';
  char diag = 'U';

  if ((p = getenv("BLAS_SIDE"))) side = *p;
  if ((p = getenv("BLAS_UPLO"))) uplo = *p;
  if ((p = getenv("BLAS_TRANS"))) trans = *p;
  if ((p = getenv("BLAS_DIAG"))) diag = *p;
  
  int m,n,k, i, j,lda,ldb;

  int from = 1;
  int to = 1;
  int step = 1;
  lc_add=0;



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

 
    if ((p = getenv("INIT_ONE"))) {
        init_to_one = 1;
    }

    if ((p = getenv("INIT_STATIC"))) {
        init_static = 1;
    }
    if ((p = getenv("BLCABC_ADD"))) {
        lc_add = 1;
    }
  
 

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

  fprintf(stderr,
          "From : %3d  To : %3d Step = %3d Side = %c Uplo = %c Trans = %c Diag "
          "= %c\n",
          from, to, step, side, uplo, trans, diag);


    if ( side=='L') {
        k = m;
    }
    else {
        k = n;
    }

    if ((a = (FLOAT *) malloc(sizeof (FLOAT) * (m + lc_add) * (k + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((b = (FLOAT *) malloc(sizeof (FLOAT) * (k + lc_add) * (n + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((bi = (FLOAT *) malloc(sizeof (FLOAT) * (k + lc_add) * (n + lc_add) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }
    if ((bref = (FLOAT *) malloc(sizeof (FLOAT) *(k + lc_add) * (n + lc_add) * COMPSIZE)) == NULL) {
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
                bi[i] = 1;
            } else {
                bi[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5;
            }
            b[i] = bi[i];
            bref[i] = bi[i];
        }
       

    }
    for (i = from; i <= to; i += step) {
        if (!has_param_m) {
            m = i;
        }
        if (!has_param_n) {
            n = i;
        }
 
        if (side== 'L') {
            k = m;
        }
        else {
            k = n;
        }
        if (trans == 'N') {
            lda = m + lc_add;
        } else {
            lda = k + lc_add;
        }
        if (trans == 'N') {
            ldb = k + lc_add;
        } else {
            ldb = n + lc_add;
        } 
        if (init_static) {
            FLOAT st = 0.0;
            for (i = 0; i < (m + lc_add) * (k + lc_add) * COMPSIZE; i++) {

                if (trans == 'N') {
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
            for (i = 0; i < (k + lc_add) * (n + lc_add) * COMPSIZE; i++) {
                if (trans == 'N') {
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
                bref[i]=b[i];
            }
            
        } else{
          for (j = 0; j < (k + lc_add) * (n + lc_add) * COMPSIZE; j++) {
            b[j] = bi[j];
            bref[j] = bi[j];
          }
        } 

        fprintf(stderr, " M=%4d, N=%4d, lda=%4d, ldb=%4d, alpha[0]=%f: \n", (int) m, (int) n, lda, ldb, alpha[0]);

    TRMM(&side, &uplo, &trans, &diag, &m, &n, alpha, a, &lda, b, &ldb);

#ifndef COMPLEX
    ref_trmm(&side, &uplo, &trans, &diag, &m, &n, alpha, a, &lda, bref, &ldb);
#else
    ref_ztrmm(&side, &uplo, &trans, &diag, &m, &n, alpha, a, &lda, bref, &ldb);
#endif
    compare_vals_aggregated(k, m * n, b, 1, bref, 1, STRINGIZE(TRMM));
  }

  return 0;
}

// void main(int argc, char *argv[]) __attribute__((weak, alias("MAIN__")));
