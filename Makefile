INCLUDES = -I. -I..
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O2 -Wall 
LDFLAGS    := -lm  -lpthread
LIB_BLAS   := ../libopenblas.a
UTIL       := bl_test.o
MFLAGS :=



	
%.o: %.c
	$(CC) $(CFLAGS) $(MFLAGS) $(INCLUDES) -c $< -o $@
	
all_double :: test_ddot test_drot test_daxpy test_dscal test_dasum test_dcopy test_dswap test_idamax  test_idamin
all_double test_ddot test_drot test_daxpy test_dscal test_dasum  test_dcopy test_dswap test_idamax test_idamin: 	MFLAGS = -DDOUBLE=1


all :: all_double 



test_ddot : dot.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dasum : asum.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_drot : rot.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_daxpy : axpy.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dscal : scal.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dcopy : copy.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dswap : swap.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_idamax : iamax.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_idamin : iamin.o $(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
clean ::
	@rm -f  *.o test_*

