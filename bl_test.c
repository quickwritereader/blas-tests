#include "bl_test.h"

static int nearlyEqual(FLOAT a, FLOAT b, double max_diff) {
    FLOAT absA = ABS(a);
    FLOAT absB = ABS(b);
    FLOAT diff = ABS(a - b);

    if (diff <= max_diff)
        return 1;

    FLOAT largest = MAX(absA, absB);

    if (diff <= (double) largest * max_diff)
        return 1;


    return 0;
}

static double getAggregatedMaxdiff(BLASLONG aggregated_N) {
    aggregated_N = MAX(1, aggregated_N);

    return EPL * (double) (2 * aggregated_N);
}

void compare_vals_zc(BLASLONG aggregated_N, BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;

    BLASLONG ind = -1;
    inc_x *= 2;
    inc_y *= 2;

    double maxDiff = getAggregatedMaxdiff(aggregated_N);


    for (; i < n; i++) {
        if (!nearlyEqual(x[ix], y[iy], maxDiff) || !nearlyEqual(x[ix + 1], y[iy + 1], maxDiff)) {
            ind = i;
            break;
        }

        ix += inc_x;
        iy += inc_y;
    }
    if (ind >= 0) {
        ERROR_LOG("%s Error Index: [%ld ] maxdiff %le \n", str, (long int) ind, maxDiff);

        int last = MIN(ind + 16, n);
        int begin = MAX(last - 16, 0);
        for (i = begin; i < last; i++)
            LOG("%ld)): {%f,%f}   {%f,%f} \n", i, x[i * inc_x ], x[i * inc_x + 1], y[i * inc_y ], y[i * inc_y + 1]);


    } else {

        PASS_LOG("%s PASSED for maxdiff %le\n", str, maxDiff);
    }
}

void compare_vals_ds(BLASLONG aggregated_N, BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;
    BLASLONG ind = -1;

    double maxDiff = getAggregatedMaxdiff(aggregated_N);

    
    for (; i < n; i++) {
        if (!nearlyEqual(x[ix], y[iy], maxDiff)) {
            ind = i;
             LOG("ERROR %ld)): {%f}   {%f } \n", i, x[ix], y[iy ]);
        }
        ix += inc_x;
        iy += inc_y;
    }
    if (ind >= 0) {
        ERROR_LOG("%s Error Index: [%ld ] maxdiff %le \n", str, (long int) ind, maxDiff);

        int last = MIN(ind + 128, n);
        int begin = MAX(last - 128, 0);
        for (i = begin; i < last; i++)
            LOG("%ld)): {%f}   {%f } \n", i, x[i * inc_x ], y[i * inc_y ]);

    } else {
        int last = MIN(ind + 128, n);
        int begin = MAX(last - 128, 0);
        for (i = begin; i < last; i++)
            LOG("%ld)): {%f}   {%f } \n", i, x[i * inc_x ], y[i * inc_y ]);
        PASS_LOG("%s PASSED for maxdiff %le\n", str, maxDiff);
    }
}

void compare_vals_aggregated(BLASLONG aggregated_N, BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str) {
#ifdef COMPLEX
    compare_vals_zc(aggregated_N, n, x, inc_x, y, inc_y, str);
#else    

    compare_vals_ds(aggregated_N, n, x, inc_x, y, inc_y, str);
#endif      
}

void compare_vals(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str) {

    compare_vals_aggregated(1, n, x, inc_x, y, inc_y, str);
}

void compare_aggregate_real(BLASLONG aggregated_N, FLOAT x, FLOAT y, const char * str) {


    double maxDiff = EPL * 8.0;
    aggregated_N = MAX(1.0, aggregated_N);

    maxDiff = EPL * (double) (2 * aggregated_N);

    if (!nearlyEqual(x, y, maxDiff)) {
        ERROR_LOG(" {%f} {%f } maxdiff %le\n", x, y, maxDiff);
    } else {
        PASS_LOG("%s PASSED for maxdiff %le\n", str, maxDiff);
    }

}

void compare_aggregate(BLASLONG aggregated_N, FLOAT *x, FLOAT *y, const char * str) {


    double maxDiff = EPL * 8.0;
    aggregated_N = MAX(1.0, aggregated_N);
    //  if (strncmp(str, "ASUM", MIN(4, strlen(str))) == 0 || strncmp(str, "DOT", MIN(3, strlen(str))) == 0) {
#ifdef COMPLEX
    maxDiff = EPL * (double) (8 * aggregated_N);
#else
    maxDiff = EPL * (double) (2 * aggregated_N);
#endif                
    //   }


#ifdef COMPLEX    
    if (!nearlyEqual(x[0], y[0], maxDiff) || !nearlyEqual(x[1], y[1], maxDiff)) {
        ERROR_LOG("{%f,%f}   {%f,%f} maxdiff %le\n", x[0 ], x[1], y[0], y[1], maxDiff);
    }
#else
    if (!nearlyEqual(x[0], y[0], maxDiff)) {
        ERROR_LOG(" {%f} {%f } maxdiff %le\n", x[0 ], y[0], maxDiff);
    }
#endif  
    else {
        PASS_LOG("%s PASSED for maxdiff %le\n", str, maxDiff);
    }

}

FLOAT ref_dot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;

    FLOAT dot = 0.0;

    if (n <= 0) return (dot);

    while (i < n) {

        dot += y[iy] * x[ix];
        ix += inc_x;
        iy += inc_y;
        i++;

    }
    return (dot);
}

COMPLEX_FLOAT ref_zcdot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y) {

    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;
    FLOAT dot[2];
    COMPLEX_FLOAT result;
    BLASLONG inc_x2;
    BLASLONG inc_y2;

    dot[0] = 0.0;
    dot[1] = 0.0;

    result.real = 0.0;
    result.imag = 0.0;

    if (n < 1) return (result);

    inc_x2 = 2 * inc_x;
    inc_y2 = 2 * inc_y;

    while (i < n) {
#if !defined(CONJ)
        dot[0] += (x[ix] * y[iy] - x[ix + 1] * y[iy + 1]);
        dot[1] += (x[ix + 1] * y[iy] + x[ix] * y[iy + 1]);
#else
        dot[0] += (x[ix] * y[iy] + x[ix + 1] * y[iy + 1]);
        dot[1] -= (x[ix + 1] * y[iy] - x[ix] * y[iy + 1]);
#endif
        ix += inc_x2;
        iy += inc_y2;
        i++;

    }
    result.real = dot[0];
    result.imag = dot[1];

    return (result);
}

FLOAT ref_asum(BLASLONG n, FLOAT *x, BLASLONG inc_x) {
    BLASLONG i = 0;
    BLASLONG ix = 0;

    FLOAT asum = 0.0;

    if (n <= 0) return (asum);
#ifdef COMPLEX
    inc_x = 2 * inc_x;

    while (i < n) {
        asum += ABS(x[ix]) + ABS(x[ix + 1]);
        ix += inc_x;
        i++;
    }
#else    
    while (i < n) {

        asum += ABS(x[ix]);
        ix += inc_x;
        i++;
    }
#endif    
    return (asum);
}

int ref_rot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;
#ifdef COMPLEX
    FLOAT temp[2];


    if (n <= 0) return (0);

    inc_x = 2 * inc_x;
    inc_y = 2 * inc_y;

    while (i < n) {
        temp[0] = c * x[ix] + s * y[iy];
        temp[1] = c * x[ix + 1] + s * y[iy + 1];
        y[iy] = c * y[iy] - s * x[ix];
        y[iy + 1] = c * y[iy + 1] - s * x[ix + 1];
        x[ix] = temp[0];
        x[ix + 1] = temp[1];

        ix += inc_x;
        iy += inc_y;
        i++;

    }
#else    
    FLOAT temp;

    if (n <= 0) return (0);


    while (i < n) {

        temp = c * x[ix] + s * y[iy];
        y[iy] = c * y[iy] - s * x[ix];
        x[ix] = temp;

        ix += inc_x;
        iy += inc_y;
        i++;

    }

#endif  
    return (0);

}

void ref_axpy(BLASLONG n, FLOAT *alpha, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;
    FLOAT da_r = alpha[0];

#ifdef COMPLEX
    FLOAT da_i = alpha[1];
    inc_x = 2 * inc_x;
    inc_y = 2 * inc_y;

    while (i < n) {
#if !defined(CONJ)
        y[iy] += (da_r * x[ix] - da_i * x[ix + 1]);
        y[iy + 1] += (da_r * x[ix + 1] + da_i * x[ix]);
#else
        y[iy] += (da_r * x[ix] + da_i * x[ix + 1]);
        y[iy + 1] -= (da_r * x[ix + 1] - da_i * x[ix]);
#endif
        ix += inc_x;
        iy += inc_y;
        i++;

    }

#else  

    while (i < n) {

        y[iy] += da_r * x[ix];
        ix += inc_x;
        iy += inc_y;
        i++;

    }

#endif    
}

void ref_scal(BLASLONG n, FLOAT *alpha, FLOAT *x, BLASLONG inc_x) {
    BLASLONG i = 0;
    BLASLONG ix = 0;
    FLOAT da_r = alpha[0];
#ifdef COMPLEX
    FLOAT temp0 = 0;
    FLOAT da_i = alpha[1];
    inc_x = 2 * inc_x;
    while (i < n) {

        temp0 = da_r * x[ix] - da_i * x[ix + 1];
        x[ix + 1] = da_r * x[ix + 1] + da_i * x[ix];
        x[ix] = temp0;
        ix += inc_x;
        i++;

    }
#else  
    while (i < n) {

        x[ix] = da_r * x[ix];
        ix += inc_x;
        i++;

    }
#endif  
}

unsigned char lsame_(char *ca, char *cb) {


    /* System generated locals */
    unsigned char ret_val;

    /* Local variables */
    BLASLONG inta, intb, zcode;


    ret_val = *(unsigned char *) ca == *(unsigned char *) cb;
    if (ret_val) {
        return ret_val;
    }

    /*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

    /*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
           machines, on which ICHAR returns a value with bit 8 set.   
           ICHAR('A') on Prime machines returns 193 which is the same as   
           ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *) ca;
    intb = *(unsigned char *) cb;

    if (zcode == 90 || zcode == 122) {

        /*        ASCII is assumed - ZCODE is the ASCII code of either lower o
        r   
                  upper case 'Z'. */

        if (inta >= 97 && inta <= 122) {
            inta += -32;
        }
        if (intb >= 97 && intb <= 122) {
            intb += -32;
        }

    } else if (zcode == 233 || zcode == 169) {

        /*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
         or   
                  upper case 'Z'. */

        if ((inta >= 129 && inta <= 137) || (inta >= 145 && inta <= 153) || (inta
                >= 162 && inta <= 169)) {
            inta += 64;
        }
        if ((intb >= 129 && intb <= 137) || (intb >= 145 && intb <= 153) || (intb
                >= 162 && intb <= 169)) {
            intb += 64;
        }

    } else if (zcode == 218 || zcode == 250) {

        /*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
        e   
                  plus 128 of either lower or upper case 'Z'. */

        if (inta >= 225 && inta <= 250) {
            inta += -32;
        }
        if (intb >= 225 && intb <= 250) {
            intb += -32;
        }
    }
    ret_val = inta == intb;

    /*     RETURN   

           End of LSAME */

    return ret_val;
} /* lsame_ */

//ref gem

int ref_gemm(char *transa, char *transb, int *m, int *
        n, int *k, FLOAT *alpha, FLOAT *a, int *lda, FLOAT *b, int *
        ldb, FLOAT *beta, FLOAT *c, int *ldc) {

    BLASLONG i__1, i__2,
            i__3;

    /* Local variables */
    BLASLONG info;
    unsigned char nota, notb;
    FLOAT temp;
    BLASLONG i, j, l, ncola;
    BLASLONG nrowa, nrowb;


    nota = lsame_(transa, "N");
    notb = lsame_(transb, "N");
    if (nota) {
        nrowa = *m;
        ncola = *k;
    } else {
        nrowa = *k;
        ncola = *m;
    }
    if (notb) {
        nrowb = *k;
    } else {
        nrowb = *n;
    }


    info = 0;
    if (!nota && !lsame_(transa, "C") && !lsame_(transa, "T")) {
        info = 1;
    } else if (!notb && !lsame_(transb, "C") && !lsame_(transb,
            "T")) {
        info = 2;
    } else if (*m < 0) {
        info = 3;
    } else if (*n < 0) {
        info = 4;
    } else if (*k < 0) {
        info = 5;
    } else if (*lda < MAX(1, nrowa)) {
        info = 8;
    } else if (*ldb < MAX(1, nrowb)) {
        info = 10;
    } else if (*ldc < MAX(1, *m)) {
        info = 13;
    }
    if (info != 0) {
        return 0;
    }

    /*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0.f || *k == 0) && *beta == 1.f) {
        return 0;
    }

    /*     And if  alpha.eq.zero. */

    if (*alpha == 0.f) {
        if (*beta == 0.f) {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    C(i, j) = 0.f;
                    /* L10: */
                }
                /* L20: */
            }
        } else {
            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    C(i, j) = *beta * C(i, j);
                    /* L30: */
                }
                /* L40: */
            }
        }
        return 0;
    }

    /*     Start the operations. */

    if (notb) {
        if (nota) {

            /*           Form  C := alpha*A*B + beta*C. */

            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (*beta == 0.f) {
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        C(i, j) = 0.f;
                        /* L50: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        C(i, j) = *beta * C(i, j);
                        /* L60: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= *k; ++l) {
                    if (B(l, j) != 0.f) {
                        temp = *alpha * B(l, j);
                        i__3 = *m;
                        for (i = 1; i <= *m; ++i) {
                            C(i, j) += temp * A(i, l);
                            /* L70: */
                        }
                    }
                    /* L80: */
                }
                /* L90: */
            }
        } else {

            /*           Form  C := alpha*A'*B + beta*C */

            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= *k; ++l) {
                        temp += A(l, i) * B(l, j);
                        /* L100: */
                    }
                    if (*beta == 0.f) {
                        C(i, j) = *alpha * temp;
                    } else {
                        C(i, j) = *alpha * temp + *beta * C(i, j);
                    }
                    /* L110: */
                }
                /* L120: */
            }
        }
    } else {
        if (nota) {

            /*           Form  C := alpha*A*B' + beta*C */

            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                if (*beta == 0.f) {
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        C(i, j) = 0.f;
                        /* L130: */
                    }
                } else if (*beta != 1.f) {
                    i__2 = *m;
                    for (i = 1; i <= *m; ++i) {
                        C(i, j) = *beta * C(i, j);
                        /* L140: */
                    }
                }
                i__2 = *k;
                for (l = 1; l <= *k; ++l) {
                    if (B(j, l) != 0.f) {
                        temp = *alpha * B(j, l);
                        i__3 = *m;
                        for (i = 1; i <= *m; ++i) {
                            C(i, j) += temp * A(i, l);
                            /* L150: */
                        }
                    }
                    /* L160: */
                }
                /* L170: */
            }
        } else {

            /*           Form  C := alpha*A'*B' + beta*C */

            i__1 = *n;
            for (j = 1; j <= *n; ++j) {
                i__2 = *m;
                for (i = 1; i <= *m; ++i) {
                    temp = 0.f;
                    i__3 = *k;
                    for (l = 1; l <= *k; ++l) {
                        temp += A(l, i) * B(j, l);
                        /* L180: */
                    }
                    if (*beta == 0.f) {
                        C(i, j) = *alpha * temp;
                    } else {
                        C(i, j) = *alpha * temp + *beta * C(i, j);
                    }
                    /* L190: */
                }
                /* L200: */
            }
        }
    }

    return 0;

    /*     End of  . */

}


void d_cnjg(COMPLEX_FLOAT *r, COMPLEX_FLOAT *z) 
{
	FLOAT zi = z->imag;
	r->real = z->real;
	r->imag = -zi; 
}

int ref_zgemm(char *transa, char *transb, int *m, int *
	n, int *k, COMPLEX_FLOAT *alpha, COMPLEX_FLOAT *a, int *lda, 
	COMPLEX_FLOAT *b, int *ldb, COMPLEX_FLOAT *beta, COMPLEX_FLOAT *c,
	 int *ldc)
{


    /* System generated locals */
    BLASLONG a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6;
    COMPLEX_FLOAT z__1, z__2, z__3, z__4;
 
 
    /* Local variables */
    BLASLONG info;
    unsigned char nota, notb,conja, conjb;
    COMPLEX_FLOAT temp;
    BLASLONG i, j, l, ncola;
    BLASLONG nrowa, nrowb;

    nota = lsame_(transa, "N");
    notb = lsame_(transb, "N");
    conja = lsame_(transa, "C");
    conjb = lsame_(transb, "C");
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }

/*     Test the input parameters. */

    info = 0;
    if (! nota && ! conja && ! lsame_(transa, "T")) {
	info = 1;
    } else if (! notb && ! conjb && ! lsame_(transb, "T")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < MAX(1,nrowa)) {
	info = 8;
    } else if (*ldb < MAX(1,nrowb)) {
	info = 10;
    } else if (*ldc < MAX(1,*m)) {
	info = 13;
    }
    if (info != 0) { 
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (alpha->real == 0. && alpha->imag == 0. || *k == 0) &&
	     (beta->real == 1. && beta->imag == 0.)) {
	return 0;
    }

/*     And when  alpha.eq.zero. */

    if (alpha->real == 0. && alpha->imag == 0.) {
	if (beta->real == 0. && beta->imag == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * c_dim1;
		    C(i,j).real = 0., C(i,j).imag = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    i__3 = i + j * c_dim1;
		    i__4 = i + j * c_dim1;
		    z__1.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
			    z__1.imag = beta->real * C(i,j).imag + beta->imag * C(i,j)
			    .real;
		    C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }

/*     Start the operations. */

    if (notb) {
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (beta->real == 0. && beta->imag == 0.) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * c_dim1;
			C(i,j).real = 0., C(i,j).imag = 0.;
/* L50: */
		    }
		} else if (beta->real != 1. || beta->imag != 0.) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			z__1.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__1.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    i__3 = l + j * b_dim1;
		    if (B(l,j).real != 0. || B(l,j).imag != 0.) {
			i__3 = l + j * b_dim1;
			z__1.real = alpha->real * B(l,j).real - alpha->imag * B(l,j).imag, 
				z__1.imag = alpha->real * B(l,j).imag + alpha->imag * B(l,j).real;
			temp.real = z__1.real, temp.imag = z__1.imag;
			i__3 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__4 = i + j * c_dim1;
			    i__5 = i + j * c_dim1;
			    i__6 = i + l * a_dim1;
			    z__2.real = temp.real * A(i,l).real - temp.imag * A(i,l).imag, 
				    z__2.imag = temp.real * A(i,l).imag + temp.imag * A(i,l).real;
			    z__1.real = C(i,j).real + z__2.real, z__1.imag = C(i,j).imag + 
				    z__2.imag;
			    C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else if (conja) {

/*           Form  C := alpha*conjg( A' )*B + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp.real = 0., temp.imag = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			d_cnjg(&z__3, &A(l,i));
			i__4 = l + j * b_dim1;
			z__2.real = z__3.real * B(l,j).real - z__3.imag * B(l,j).imag, 
				z__2.imag = z__3.real * B(l,j).imag + z__3.imag * B(l,j)
				.real;
			z__1.real = temp.real + z__2.real, z__1.imag = temp.imag + z__2.imag;
			temp.real = z__1.real, temp.imag = z__1.imag;
/* L100: */
		    }
		    if (beta->real == 0. && beta->imag == 0.) {
			i__3 = i + j * c_dim1;
			z__1.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__1.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__2.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			i__4 = i + j * c_dim1;
			z__3.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__3.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    }
/* L110: */
		}
/* L120: */
	    }
	} else {

/*           Form  C := alpha*A'*B + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp.real = 0., temp.imag = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			i__4 = l + i * a_dim1;
			i__5 = l + j * b_dim1;
			z__2.real = A(l,i).real * B(l,j).real - A(l,i).imag * B(l,j)
				.imag, z__2.imag = A(l,i).real * B(l,j).imag + A(l,i)
				.imag * B(l,j).real;
			z__1.real = temp.real + z__2.real, z__1.imag = temp.imag + z__2.imag;
			temp.real = z__1.real, temp.imag = z__1.imag;
/* L130: */
		    }
		    if (beta->real == 0. && beta->imag == 0.) {
			i__3 = i + j * c_dim1;
			z__1.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__1.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__2.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			i__4 = i + j * c_dim1;
			z__3.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__3.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    }
/* L140: */
		}
/* L150: */
	    }
	}
    } else if (nota) {
	if (conjb) {

/*           Form  C := alpha*A*conjg( B' ) + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (beta->real == 0. && beta->imag == 0.) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * c_dim1;
			C(i,j).real = 0., C(i,j).imag = 0.;
/* L160: */
		    }
		} else if (beta->real != 1. || beta->imag != 0.) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			z__1.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__1.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
/* L170: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    i__3 = j + l * b_dim1;
		    if (B(j,l).real != 0. || B(j,l).imag != 0.) {
			d_cnjg(&z__2, &B(j,l));
			z__1.real = alpha->real * z__2.real - alpha->imag * z__2.imag, 
				z__1.imag = alpha->real * z__2.imag + alpha->imag * 
				z__2.real;
			temp.real = z__1.real, temp.imag = z__1.imag;
			i__3 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__4 = i + j * c_dim1;
			    i__5 = i + j * c_dim1;
			    i__6 = i + l * a_dim1;
			    z__2.real = temp.real * A(i,l).real - temp.imag * A(i,l).imag, 
				    z__2.imag = temp.real * A(i,l).imag + temp.imag * A(i,l).real;
			    z__1.real = C(i,j).real + z__2.real, z__1.imag = C(i,j).imag + 
				    z__2.imag;
			    C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
/* L180: */
			}
		    }
/* L190: */
		}
/* L200: */
	    }
	} else {

/*           Form  C := alpha*A*B'          + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (beta->real == 0. && beta->imag == 0.) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * c_dim1;
			C(i,j).real = 0., C(i,j).imag = 0.;
/* L210: */
		    }
		} else if (beta->real != 1. || beta->imag != 0.) {
		    i__2 = *m;
		    for (i = 1; i <= *m; ++i) {
			i__3 = i + j * c_dim1;
			i__4 = i + j * c_dim1;
			z__1.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__1.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
/* L220: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= *k; ++l) {
		    i__3 = j + l * b_dim1;
		    if (B(j,l).real != 0. || B(j,l).imag != 0.) {
			i__3 = j + l * b_dim1;
			z__1.real = alpha->real * B(j,l).real - alpha->imag * B(j,l).imag, 
				z__1.imag = alpha->real * B(j,l).imag + alpha->imag * B(j,l).real;
			temp.real = z__1.real, temp.imag = z__1.imag;
			i__3 = *m;
			for (i = 1; i <= *m; ++i) {
			    i__4 = i + j * c_dim1;
			    i__5 = i + j * c_dim1;
			    i__6 = i + l * a_dim1;
			    z__2.real = temp.real * A(i,l).real - temp.imag * A(i,l).imag, 
				    z__2.imag = temp.real * A(i,l).imag + temp.imag * A(i,l).real;
			    z__1.real = C(i,j).real + z__2.real, z__1.imag = C(i,j).imag + 
				    z__2.imag;
			    C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
/* L230: */
			}
		    }
/* L240: */
		}
/* L250: */
	    }
	}
    } else if (conja) {
	if (conjb) {

/*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp.real = 0., temp.imag = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			d_cnjg(&z__3, &A(l,i));
			d_cnjg(&z__4, &B(j,l));
			z__2.real = z__3.real * z__4.real - z__3.imag * z__4.imag, z__2.imag = 
				z__3.real * z__4.imag + z__3.imag * z__4.real;
			z__1.real = temp.real + z__2.real, z__1.imag = temp.imag + z__2.imag;
			temp.real = z__1.real, temp.imag = z__1.imag;
/* L260: */
		    }
		    if (beta->real == 0. && beta->imag == 0.) {
			i__3 = i + j * c_dim1;
			z__1.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__1.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__2.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			i__4 = i + j * c_dim1;
			z__3.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__3.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    }
/* L270: */
		}
/* L280: */
	    }
	} else {

/*           Form  C := alpha*conjg( A' )*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp.real = 0., temp.imag = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			d_cnjg(&z__3, &A(l,i));
			i__4 = j + l * b_dim1;
			z__2.real = z__3.real * B(j,l).real - z__3.imag * B(j,l).imag, 
				z__2.imag = z__3.real * B(j,l).imag + z__3.imag * B(j,l)
				.real;
			z__1.real = temp.real + z__2.real, z__1.imag = temp.imag + z__2.imag;
			temp.real = z__1.real, temp.imag = z__1.imag;
/* L290: */
		    }
		    if (beta->real == 0. && beta->imag == 0.) {
			i__3 = i + j * c_dim1;
			z__1.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__1.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__2.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			i__4 = i + j * c_dim1;
			z__3.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__3.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    }
/* L300: */
		}
/* L310: */
	    }
	}
    } else {
	if (conjb) {

/*           Form  C := alpha*A'*conjg( B' ) + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp.real = 0., temp.imag = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			i__4 = l + i * a_dim1;
			d_cnjg(&z__3, &B(j,l));
			z__2.real = A(l,i).real * z__3.real - A(l,i).imag * z__3.imag, 
				z__2.imag = A(l,i).real * z__3.imag + A(l,i).imag * 
				z__3.real;
			z__1.real = temp.real + z__2.real, z__1.imag = temp.imag + z__2.imag;
			temp.real = z__1.real, temp.imag = z__1.imag;
/* L320: */
		    }
		    if (beta->real == 0. && beta->imag == 0.) {
			i__3 = i + j * c_dim1;
			z__1.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__1.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__2.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			i__4 = i + j * c_dim1;
			z__3.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__3.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    }
/* L330: */
		}
/* L340: */
	    }
	} else {

/*           Form  C := alpha*A'*B' + beta*C */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    temp.real = 0., temp.imag = 0.;
		    i__3 = *k;
		    for (l = 1; l <= *k; ++l) {
			i__4 = l + i * a_dim1;
			i__5 = j + l * b_dim1;
			z__2.real = A(l,i).real * B(j,l).real - A(l,i).imag * B(j,l)
				.imag, z__2.imag = A(l,i).real * B(j,l).imag + A(l,i)
				.imag * B(j,l).real;
			z__1.real = temp.real + z__2.real, z__1.imag = temp.imag + z__2.imag;
			temp.real = z__1.real, temp.imag = z__1.imag;
/* L350: */
		    }
		    if (beta->real == 0. && beta->imag == 0.) {
			i__3 = i + j * c_dim1;
			z__1.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__1.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    } else {
			i__3 = i + j * c_dim1;
			z__2.real = alpha->real * temp.real - alpha->imag * temp.imag, 
				z__2.imag = alpha->real * temp.imag + alpha->imag * 
				temp.real;
			i__4 = i + j * c_dim1;
			z__3.real = beta->real * C(i,j).real - beta->imag * C(i,j).imag, 
				z__3.imag = beta->real * C(i,j).imag + beta->imag * C(i,j).real;
			z__1.real = z__2.real + z__3.real, z__1.imag = z__2.imag + z__3.imag;
			C(i,j).real = z__1.real, C(i,j).imag = z__1.imag;
		    }
/* L360: */
		}
/* L370: */
	    }
	}
    }

    return 0;

/*     End of ZGEMM . */

} /* zgemm_ */



int ref_trmm(char* side, char* uplo, char* transa, char* diag,
    int* m, int* n, FLOAT* alpha, FLOAT* a, int* lda, FLOAT* b,
    int* ldb)
{
    /* System generated locals */
    BLASLONG a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    BLASLONG i__, j, k, info;
    FLOAT temp;
    unsigned char  lside;
    BLASLONG nrowa;
    unsigned char  upper;

    unsigned char  nounit;

    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    /* Function Body */
    lside = lsame_(side, "L");
    if (lside) {
        nrowa = *m;
    }
    else {
        nrowa = *n;
    }
    nounit = lsame_(diag, "N");
    upper = lsame_(uplo, "U");
    info = 0;
    if (!lside && !lsame_(side, "R")) {
        info = 1;
    }
    else if (!upper && !lsame_(uplo, "L")) {
        info = 2;
    }
    else if (!lsame_(transa, "N") && !lsame_(transa, "T") && !lsame_(transa, "C")) {
        info = 3;
    }
    else if (!lsame_(diag, "U") && !lsame_(diag, "N")) {
        info = 4;
    }
    else if (*m < 0) {
        info = 5;
    }
    else if (*n < 0) {
        info = 6;
    }
    else if (*lda < MAX(1, nrowa)) {
        info = 9;
    }
    else if (*ldb < MAX(1, *m)) {
        info = 11;
    }
    if (info != 0) {

        return 0;
    }
    /*     Quick return if possible. */
    if (*n == 0) {
        return 0;
    }
    /*     And when  alpha.eq.zero. */
    if (*alpha == 0.f) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *m;
            for (i__ = 1; i__ <= i__2; ++i__) {
                b[i__ + j * b_dim1] = 0.f;
                /* L10: */
            }
            /* L20: */
        }
        return 0;
    }
    /*     Start the operations. */
    if (lside) {
        if (lsame_(transa, "N")) {
            /*           Form  B := alpha*A*B. */
            if (upper) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (k = 1; k <= i__2; ++k) {
                        if (b[k + j * b_dim1] != 0.f) {
                            temp = *alpha * b[k + j * b_dim1];
                            i__3 = k - 1;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * a[i__ + k * a_dim1];
                                /* L30: */
                            }
                            if (nounit) {
                                temp *= a[k + k * a_dim1];
                            }
                            b[k + j * b_dim1] = temp;
                        }
                        /* L40: */
                    }
                    /* L50: */
                }
            }
            else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    for (k = *m; k >= 1; --k) {
                        if (b[k + j * b_dim1] != 0.f) {
                            temp = *alpha * b[k + j * b_dim1];
                            b[k + j * b_dim1] = temp;
                            if (nounit) {
                                b[k + j * b_dim1] *= a[k + k * a_dim1];
                            }
                            i__2 = *m;
                            for (i__ = k + 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * a[i__ + k * a_dim1];
                                /* L60: */
                            }
                        }
                        /* L70: */
                    }
                    /* L80: */
                }
            }
        }
        else {
            /*           Form  B := alpha*A'*B. */
            if (upper) {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    for (i__ = *m; i__ >= 1; --i__) {
                        temp = b[i__ + j * b_dim1];
                        if (nounit) {
                            temp *= a[i__ + i__ * a_dim1];
                        }
                        i__2 = i__ - 1;
                        for (k = 1; k <= i__2; ++k) {
                            temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
                            /* L90: */
                        }
                        b[i__ + j * b_dim1] = *alpha * temp;
                        /* L100: */
                    }
                    /* L110: */
                }
            }
            else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        temp = b[i__ + j * b_dim1];
                        if (nounit) {
                            temp *= a[i__ + i__ * a_dim1];
                        }
                        i__3 = *m;
                        for (k = i__ + 1; k <= i__3; ++k) {
                            temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
                            /* L120: */
                        }
                        b[i__ + j * b_dim1] = *alpha * temp;
                        /* L130: */
                    }
                    /* L140: */
                }
            }
        }
    }
    else {
        if (lsame_(transa, "N")) {
            /*           Form  B := alpha*B*A. */
            if (upper) {
                for (j = *n; j >= 1; --j) {
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__1 = *m;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        /* L150: */
                    }
                    i__1 = j - 1;
                    for (k = 1; k <= i__1; ++k) {
                        if (a[k + j * a_dim1] != 0.f) {
                            temp = *alpha * a[k + j * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
                                /* L160: */
                            }
                        }
                        /* L170: */
                    }
                    /* L180: */
                }
            }
            else {
                i__1 = *n;
                for (j = 1; j <= i__1; ++j) {
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[j + j * a_dim1];
                    }
                    i__2 = *m;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
                        /* L190: */
                    }
                    i__2 = *n;
                    for (k = j + 1; k <= i__2; ++k) {
                        if (a[k + j * a_dim1] != 0.f) {
                            temp = *alpha * a[k + j * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
                                /* L200: */
                            }
                        }
                        /* L210: */
                    }
                    /* L220: */
                }
            }
        }
        else {
            /*           Form  B := alpha*B*A'. */
            if (upper) {
                i__1 = *n;
                for (k = 1; k <= i__1; ++k) {
                    i__2 = k - 1;
                    for (j = 1; j <= i__2; ++j) {
                        if (a[j + k * a_dim1] != 0.f) {
                            temp = *alpha * a[j + k * a_dim1];
                            i__3 = *m;
                            for (i__ = 1; i__ <= i__3; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
                                /* L230: */
                            }
                        }
                        /* L240: */
                    }
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[k + k * a_dim1];
                    }
                    if (temp != 1.f) {
                        i__2 = *m;
                        for (i__ = 1; i__ <= i__2; ++i__) {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                            /* L250: */
                        }
                    }
                    /* L260: */
                }
            }
            else {
                for (k = *n; k >= 1; --k) {
                    i__1 = *n;
                    for (j = k + 1; j <= i__1; ++j) {
                        if (a[j + k * a_dim1] != 0.f) {
                            temp = *alpha * a[j + k * a_dim1];
                            i__2 = *m;
                            for (i__ = 1; i__ <= i__2; ++i__) {
                                b[i__ + j * b_dim1] += temp * b[i__ + k * b_dim1];
                                /* L270: */
                            }
                        }
                        /* L280: */
                    }
                    temp = *alpha;
                    if (nounit) {
                        temp *= a[k + k * a_dim1];
                    }
                    if (temp != 1.f) {
                        i__1 = *m;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
                            /* L290: */
                        }
                    }
                    /* L300: */
                }
            }
        }
    }
    return 0;
    /*     End of ref_TRMM . */
} /* ref_trmm_ */
 