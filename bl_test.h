#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <cblas.h>
#ifndef BL_TEST_H
#define BL_TEST_H



#ifdef __cplusplus
extern "C" {
#endif





#ifndef MIN

#define MIN(a,b)   (a>b? b:a)

#endif


#define MAX(a,b)   (a<b? b:a)



#ifdef XDOUBLE
#define FLOAT xdouble
#ifdef QUAD_PRECISION
#define XFLOAT xidouble
#endif
#ifdef QUAD_PRECISION
#define SIZE 32
#define  BASE_SHIFT 5
#define ZBASE_SHIFT 6
#else
#define SIZE 16
#define  BASE_SHIFT 4
#define ZBASE_SHIFT 5
#endif
#elif defined(DOUBLE)
#define FLOAT double
#define SIZE 8
#define  BASE_SHIFT 3
#define ZBASE_SHIFT 4
#else
#define FLOAT float
#define SIZE    4
#define  BASE_SHIFT 2
#define ZBASE_SHIFT 3
#endif

#ifndef XFLOAT
#define XFLOAT FLOAT
#endif

#ifndef COMPLEX
#define COMPSIZE  1
 
#else
#define COMPSIZE  2    
#endif

#ifdef DOUBLE 
#define ABS fabs
#define FLOAT_MIN_NORMAL DBL_MIN
#define  FLOAT_MAX_VALUE DBL_MAX 
#define EPL DBL_EPSILON 
#else 
#define ABS fabsf
#define FLOAT_MIN_NORMAL FLT_MIN
#define  FLOAT_MAX_VALUE FLT_MAX 
#define EPL FLT_EPSILON     
#endif
#define CABS1(x,i)    ABS(x[i])+ABS(x[i+1])

#define STRINGIZE(x) #x
#define LOG(format, ...)  fprintf (stderr, format  , ##__VA_ARGS__)
#define ERROR_LOG(format, ...)  fprintf (stderr,"Error: " format "\n", ##__VA_ARGS__)
#define PASS_LOG(format, ...)  fprintf (stderr,"Passed: " format "\n", ##__VA_ARGS__)    
    typedef struct __attribute__((__packed__)) {
        FLOAT real, imag;
    } COMPLEX_FLOAT;

    void compare_aggregate(BLASLONG aggregated_N,FLOAT *x,FLOAT *y, const char * str);
    void compare_aggregate_real(BLASLONG aggregated_N, FLOAT  x, FLOAT  y, const char * str);
    void compare_vals_ds(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str);
    void compare_vals_zc(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str);
        
    void compare_vals(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str);

    /*Reference functions*/
    FLOAT ref_dot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y);

    COMPLEX_FLOAT ref_zcdot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y);

    FLOAT ref_asum(BLASLONG n, FLOAT *x, BLASLONG inc_x);

    int ref_rot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s);

    void ref_axpy(BLASLONG n, FLOAT *alpha, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y);

    void ref_scal(BLASLONG n, FLOAT *alpha, FLOAT *x, BLASLONG inc_x);

#ifdef __cplusplus
}
#endif

#endif /* BL_TEST_H */

