
#include <bl_test.h>

#undef IAMIN



#ifdef COMPLEX
#ifdef DOUBLE
#define IAMIN   BLASFUNC(izamin)
#else
#define IAMIN   BLASFUNC(icamin)
#endif
#else
#ifdef DOUBLE
#define IAMIN   BLASFUNC(idamin)
#else
#define IAMIN   BLASFUNC(isamin)
#endif
#endif
int main(int argc, char *argv[]) {

    FLOAT *x  ; 
    FLOAT minf;
    blasint  i,ix,min,ret_min;
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
    
    if (argc > 0) {
        N = atol(*argv);
        argc--;
        argv++;
    }
     


    if ((p = getenv("OPENBLAS_INCX"))){
        inc_x = atoi(p);
        if(inc_x<=0) inc_x=1;
        LOG( "inc_x reset to 1\n");      
    }
 

    LOG( "N: %d Inc_x = %d  \n", N, inc_x );

    if ((x = (FLOAT *) malloc(sizeof (FLOAT) * N * abs(inc_x) * COMPSIZE)) == NULL) {
        LOG( "Out of Memory!!\n");
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

    min=0;
        //Find real Min Value and Index
#ifdef COMPLEX 
    inc_x2=2*inc_x;
    minf = CABS1(x,0);
    ix += inc_x2;
    i++;

    while(i < N)
    {
        if( CABS1(x,ix) < minf )
        {
            min = i;
            minf = CABS1(x,ix);
        }
        ix += inc_x2;
        i++;
    }
#else   
    minf = ABS(x[i]);
    ix += inc_x;
    i++;
while (i < N) {
            if (ABS(x[ix]) < minf) {
                min = i;
                minf = ABS(x[ix]);
            }
            ix += inc_x;
            i++;
        }
#endif       
        LOG( " %6d : Compare min \n", (int) N);
        int not_passed=0;
        ret_min=IAMIN (&N, x, &inc_x);
        LOG( "%d %c= %d\t",ret_min,(ret_min==min+1)?'=':'!',min+1);        
        not_passed=(ret_min==min+1)?0:1;
        LOG( "\n %6d : Compare shifting and swapping minimum index from 0 to MIN(N,127) \n", (int) N);
        
        for(i=0;i<MIN(N,127);i++){
            //swap
#ifdef COMPLEX
            FLOAT temp[2];
            temp[0]  = x[i*inc_x2]   ;
            temp[1]  = x[i*inc_x2+1] ;
            x[i*inc_x2]    = x[min*inc_x2]   ;
            x[i*inc_x2+1]  = x[min*inc_x2+1] ;
            x[min*inc_x2]    = temp[0] ;
            x[min*inc_x2+1]  = temp[1] ;
            min=i;
#else
            FLOAT temp=x[i*inc_x];
            x[i*inc_x] = x[min*inc_x];
            x[min*inc_x] = temp;
            min=i;
#endif
           ret_min=IAMIN (&N, x, &inc_x);
           not_passed=(ret_min==min+1)?0:1;
           LOG(  "%d %c= %d\t",ret_min,(ret_min==min+1)?'=':'!',min+1); 
            
            
        }
        LOG("\n");
        LOG("\n TEST with duplicate min values to check minimum index:\n NOTE: important for simd vector code case  \n");

        int minnx=MIN(N,127);

        int duplicate_begin=MAX(0,minnx-13);
printf("%d\n",duplicate_begin);
        for(i=duplicate_begin;i<minnx;i++){
#ifdef COMPLEX
          x[i*inc_x2]=x[min*inc_x2];
          x[i*inc_x2+1]=x[min*inc_x2+1];
#else
          x[i*inc_x]=x[min*inc_x];

#endif

        }

           ret_min=IAMIN (&N, x, &inc_x);
        printf("%d minnx %d duplicate_begin %d\n",N,minnx,duplicate_begin);

        not_passed |=(ret_min==duplicate_begin+1)?0:1;
           LOG( "%d %c= %d \t",ret_min,(ret_min==duplicate_begin+1)?'=':'!',duplicate_begin+1);


          i=duplicate_begin=MAX(0,minnx-14);

#ifdef COMPLEX
          x[i*inc_x2]=x[min*inc_x2];
          x[i*inc_x2+1]=x[min*inc_x2+1];
#else
          x[i*inc_x]=x[min*inc_x];

#endif
           ret_min=IAMIN (&N, x, &inc_x);
            not_passed |=(ret_min==duplicate_begin+1)?0:1;
           LOG( "%d %c= %d \t",ret_min,(ret_min==duplicate_begin+1)?'=':'!',duplicate_begin+1);


         LOG("\n");
        if(not_passed) {
          ERROR_LOG("%s",STRINGIZE(IAMIN));
        }else{
           PASS_LOG("%s",STRINGIZE(IAMIN));
            
        }
        LOG(  "------------\n");

        free(x);

return 0;
}
 

    	

    	 
