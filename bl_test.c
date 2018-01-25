#include "bl_test.h"

void compare_vals(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, const char * str){
    BLASLONG i=0;
    BLASLONG ix=0,iy=0;    
    double absolute_error=-1000.0;
    BLASLONG ind=-1;
    for(;i<n;i++){
        FLOAT err =(double) ABS(y[iy]-x[ix]);
        if(err>absolute_error){
            ind=i;
            absolute_error=err;
        }
        ix  += inc_x ;
        iy  += inc_y ;
    }
    fprintf(stderr,"%s MaxError Index: [%ld] Abs: %le  Rel: %.6f \n",str,(long int)ind,absolute_error,absolute_error/ABS(y[ind*inc_y]));
}

 
FLOAT ref_dot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y)
{
    BLASLONG i=0;
    BLASLONG ix=0,iy=0;

    FLOAT  dot = 0.0 ;

    if ( n <= 0 )  return(dot); 

    while(i < n)
    {

        dot += y[iy] * x[ix] ;
        ix  += inc_x ;
        iy  += inc_y ;
        i++ ;

    } 
    return(dot); 
}

 
FLOAT ref_asum(BLASLONG n, FLOAT *x, BLASLONG inc_x)
{
    BLASLONG i=0;
    BLASLONG ix=0;

    FLOAT  asum = 0.0 ;

    if ( n <= 0 )  return(asum); 

    while(i < n)
    {

        asum += ABS(x[ix]) ;
        ix  += inc_x ; 
        i++ ;

    } 
    return(asum); 
}

int ref_rot(BLASLONG n, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y, FLOAT c, FLOAT s)
{
    BLASLONG i=0;
    BLASLONG ix=0,iy=0;
     
    FLOAT temp;

    if ( n <= 0 )  return(0);

 
        while(i < n)
        {
            temp   = c*x[ix] + s*y[iy] ;
            y[iy]  = c*y[iy] - s*x[ix] ;
            x[ix]  = temp ;

            ix += inc_x ;
            iy += inc_y ;
            i++ ;

        }

   
    return(0);

}

void ref_axpy(BLASLONG n, FLOAT alpha, FLOAT *x, BLASLONG inc_x, FLOAT *y, BLASLONG inc_y ){
  BLASLONG i=0;
  BLASLONG ix=0,iy=0;    
  while(i < n)
    {

        y[iy] += alpha * x[ix] ;
        ix  += inc_x ;
        iy  += inc_y ;
        i++ ;

    }   
}

void ref_scal(BLASLONG n, FLOAT alpha, FLOAT *x, BLASLONG inc_x ){
  BLASLONG i=0;
  BLASLONG ix=0 ;    
  while(i < n)
    {

        x[ix] = alpha * x[ix] ;
        ix  += inc_x ; 
        i++ ;

    }   
}