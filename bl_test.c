#include "bl_test.h"

void compare_vals_zc(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;
    double absolute_error = -1000.0;
    BLASLONG ind = -1; 
    inc_x *= 2;
    inc_y *= 2;
    int is_image = 0;
    for (; i < n; i++) {
        FLOAT err = (double) ABS(y[iy] - x[ix]);
        if (err > absolute_error) {
            ind = i;
            absolute_error = err;
            is_image = 0;
        }
        err = (double) ABS(y[iy + 1] - x[ix + 1]);
        if (err > absolute_error) {
            ind = i;
            absolute_error = err;
            is_image = 1;
        }
        ix += inc_x;
        iy += inc_y;
    }
    fprintf(stderr, "%s MaxError Index: [%ld{%c}]  Abs: %le  Rel: %.6f \n", str, (long int) ind, is_image ? 'i' : 'r', absolute_error, absolute_error / ABS(y[ind * inc_y + is_image]));
 
}

void compare_vals_ds(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str) {
    BLASLONG i = 0;
    BLASLONG ix = 0, iy = 0;
    double absolute_error = -1000.0;
    BLASLONG ind = -1;
  
    for (; i < n; i++) {
        FLOAT err = (double) ABS(y[iy] - x[ix]);
        if (err > absolute_error) {
            ind = i;
            absolute_error = err;
        }
        ix += inc_x;
        iy += inc_y;
    }
    fprintf(stderr, "%s MaxError Index: [%ld] Abs: %le  Rel: %.6f \n", str, (long int) ind, absolute_error, absolute_error / ABS(y[ind * inc_y]));
   
}

void compare_vals(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str) {
 
#ifdef COMPLEX
    compare_vals_zc( n,  x,  inc_x,  y,  inc_y,   str) ;
#else    
   compare_vals_ds( n,  x,  inc_x,  y,  inc_y,   str) ;
#endif    
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