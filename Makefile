
BLAS_INC_DIR := $(if $(BLAS_INC_DIR),$(BLAS_INC_DIR),/opt/OpenBlas/include) 
BLAS_LIB_DIR :=  $(if $(BLAS_LIB_DIR),$(BLAS_LIB_DIR),/opt/OpenBlas/lib) 
BLAS_LIBS  := $(if $(BLAS_LIBS),  $(BLAS_LIBS),libopenblas.a )
INCLUDES = -I. -I$(BLAS_INC_DIR)
CC	   := gcc
LINKER   := $(CC)
CFLAGS   :=  -no-pie  -g -Wall -fopenmp
LDFLAGS := -lm  -lpthread
LIB_BLAS  :=  -L${BLAS_LIB_DIR}  $(foreach libx,$(BLAS_LIBS), -l:$(libx))
MFLAGS :=
PREF  :=

Double_Tests := test_ddot test_drot test_daxpy test_dscal test_dasum test_dcopy test_dswap test_idamax  test_idamin test_dgemm test_dtrmm
Float_Tests  :=  $(foreach tests,$(Double_Tests), $(subst _id,_is,$(subst _d,_s,$(tests)))) 
Complex_Double_Tests :=   $(foreach tests,$(Double_Tests), $(subst _id,_iz,$(subst _d,_z,$(tests)))) 
Complex_Float_Tests :=   $(foreach tests,$(Double_Tests), $(subst _id,_ic,$(subst _d,_c,$(tests)))) 

.PHONY:  clean


double :: $(Double_Tests) 
single :: $(Float_Tests)
zdouble :: $(Complex_Double_Tests) 
csingle :: $(Complex_Float_Tests)
all :: single
all :: double
all :: csingle
all :: zdouble

SOURCES=$(wildcard *.c)
OBJECTS=$(SOURCES:%.c=%.o)

d%.o: %.c
	$(CC) $(CFLAGS) -DDOUBLE=1 $(INCLUDES) -c $< -o $@
s%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@   
c%.o: %.c
	$(CC) $(CFLAGS) -DCOMLEX=1 $(INCLUDES) -c $< -o $@
z%.o: %.c
	$(CC) $(CFLAGS) -DDOUBLE=1 -DCOMPLEX=1 $(INCLUDES) -c $< -o $@	 
	
test_ddot : ddot.o  dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dasum : dasum.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_drot : drot.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_daxpy : daxpy.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dscal : dscal.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dcopy : dcopy.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dswap : dswap.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_idamax : diamax.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_idamin : diamin.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dgemm : dgemm.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dtrmm : dtrmm.o dbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)		
test_sdot : sdot.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_sasum : sasum.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_srot : srot.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_saxpy : saxpy.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_sscal : sscal.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_scopy : scopy.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_sswap : sswap.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_isamax : siamax.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_isamin : siamin.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)  
test_sgemm : sgemm.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_strmm : strmm.o sbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_zdot : zdot.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zasum : zasum.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_zrot : zrot.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_zaxpy : zaxpy.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zscal : zscal.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zcopy : zcopy.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zswap : zswap.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_izamax : ziamax.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_izamin : ziamin.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zgemm : zgemm.o zbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_ztrmm : ztrmm.o zbl_test.o
	echo "ztrmm check  is missin"
test_cdot : cdot.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_casum : casum.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_crot : crot.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_caxpy : caxpy.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_cscal : cscal.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_ccopy : ccopy.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_cswap : cswap.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_icamax : ciamax.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_icamin : ciamin.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS) 
test_cgemm : cgemm.o cbl_test.o
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_ctrmm : ctrmm.o cbl_test.o
	echo "ctrmm check is missing"
clean ::
	@rm -f  *.o test_*
