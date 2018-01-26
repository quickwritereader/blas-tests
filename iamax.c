
#include <bl_test.h>

#undef IAMAX



#ifdef COMPLEX
#ifdef DOUBLE
#define IAMAX   BLASFUNC(izamax)
#else
#define IAMAX   BLASFUNC(icamax)
#endif
#else
#ifdef DOUBLE
#define IAMAX   BLASFUNC(idamax)
#else
#define IAMAX   BLASFUNC(isamax)
#endif
#endif
int main(int argc, char *argv[]) {

    FLOAT *x  ; 
    FLOAT maxf;
    blasint  i,ix,max,ret_max;
    blasint inc_x = 1 ;
#if COMPLEX    
    blasint inc_x2 = 2 ;
#endif    
    char *p;

    int N = 1;
 
    argc--;
    argv++;

    if (argc > 0) {
        N = atol(*argv);
        argc--;
        argv++;
    }
     


    if ((p = getenv("OPENBLAS_INCX"))){
        inc_x = atoi(p);
        if(inc_x<=0) inc_x=1;
        fprintf(stderr, "inc_x reset to 1\n");      
    }
 

    fprintf(stderr, "N: %d Inc_x = %d  \n", N, inc_x );

    if ((x = (FLOAT *) malloc(sizeof (FLOAT) * N * abs(inc_x) * COMPSIZE)) == NULL) {
        fprintf(stderr, "Out of Memory!!\n");
        exit(1);
    }

 

   

#ifdef linux
    srandom(getpid());
#endif

 

    for (i = 0; i < N * COMPSIZE * abs(inc_x); i++) {
            x[i] = ((FLOAT) rand() / (FLOAT) RAND_MAX) - 0.5; 
        }
 
    i=0;
    ix=0;
    inc_x=inc_x>0?inc_x:-inc_x;

    max=0;
        //Find real Max Value and Index
#ifdef COMPLEX 
    inc_x2=2*inc_x;
    maxf = CABS1(x,0);
    ix += inc_x2;
    i++;

    while(i < N)
    {
        if( CABS1(x,ix) > maxf )
        {
            max = i;
            maxf = CABS1(x,ix);
        }
        ix += inc_x2;
        i++;
    }
#else   
    maxf = ABS(x[i]);
    ix += inc_x;
    i++;
while (i < N) {
            if (ABS(x[ix]) > maxf) {
                max = i;
                maxf = ABS(x[ix]);
            }
            ix += inc_x;
            i++;
        }
#endif       
        fprintf(stderr, " %6d : Compare max \n", (int) N);
        
        ret_max=IAMAX (&N, x, &inc_x);
        fprintf(stderr, "%d %c= %d\t",ret_max,(ret_max==max+1)?'=':'!',max+1);        
        
        fprintf(stderr, "\n %6d : Compare shifting and swapping maximum index from 0 to MIN(N,127) \n", (int) N);
        
        for(i=0;i<MIN(N,127);i++){
            //swap
#ifdef COMPLEX
            FLOAT temp[2];
            temp[0]  = x[i*inc_x2]   ;
            temp[1]  = x[i*inc_x2+1] ;
            x[i*inc_x2]    = x[max*inc_x2]   ;
            x[i*inc_x2+1]  = x[max*inc_x2+1] ;
            x[max*inc_x2]    = temp[0] ;
            x[max*inc_x2+1]  = temp[1] ;
            max=i;
#else
            FLOAT temp=x[i*inc_x];
            x[i*inc_x] = x[max*inc_x];
            x[max*inc_x] = temp;
            max=i;
#endif
           ret_max=IAMAX (&N, x, &inc_x);
           fprintf(stderr, "%d %c= %d\t",ret_max,(ret_max==max+1)?'=':'!',max+1);  
            
            
        }
         
        fprintf(stderr, "------------\n");

        free(x);

return 0;
}
 

    	

    	 