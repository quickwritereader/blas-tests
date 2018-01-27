
BLAS_INC_DIR = ..
BLAS_LIB_DIR = ..
BLAS_LIBS := libopenblas.a
INCLUDES = -I. -I$(BLAS_INC_DIR)
CC	   := gcc
LINKER   := $(CC)
CFLAGS   := -O2 -Wall 
LDFLAGS := -lm  -lpthread
LIB_BLAS   := $(foreach libx,$(BLAS_LIBS), $(BLAS_LIB_DIR)/$(libx))
UTIL	   := bl_test.o
MFLAGS :=
PREF  :=

Double_Tests := test_ddot test_drot test_daxpy test_dscal test_dasum test_dcopy test_dswap test_idamax  test_idamin
Float_Tests  :=  $(foreach tests,$(Double_Tests), $(subst _id,_is,$(subst _d,_s,$(tests)))) 
Complex_Double_Tests :=   $(foreach tests,$(Double_Tests), $(subst _id,_iz,$(subst _d,_z,$(tests)))) 
Complex_Float_Tests :=   $(foreach tests,$(Double_Tests), $(subst _id,_ic,$(subst _d,_c,$(tests)))) 

.PHONY:  clean



double :: $(Double_Tests)
double $(Double_Tests):  MFLAGS = -DDOUBLE=1 
single :: $(Float_Tests)
single $(Float_Tests) :  MFLAGS= 
zdouble :: $(Complex_Double_Tests)
zdouble $(Complex_Double_Tests):	 MFLAGS =  -DCOMPLEX=1 -DDOUBLE=1  
csingle :: $(Complex_Float_Tests)
csingle $(Complex_Float_Tests) :	 MFLAGS=  -DCOMPLEX=1
all :: single
all :: double
all :: csingle
all :: zdouble

SOURCES=$(wildcard *.c)
OBJECTS=$(SOURCES:%.c=%.o)

d%.o: %.c
	$(CC) $(CFLAGS) $(MFLAGS) $(INCLUDES) -c $< -o $@
s%.o: %.c
	$(CC) $(CFLAGS) $(MFLAGS) $(INCLUDES) -c $< -o $@   
c%.o: %.c
	$(CC) $(CFLAGS) $(MFLAGS) $(INCLUDES) -c $< -o $@
z%.o: %.c
	$(CC) $(CFLAGS) $(MFLAGS) $(INCLUDES) -c $< -o $@	 
	

test_ddot : ddot.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dasum : dasum.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_drot : drot.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_daxpy : daxpy.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dscal : dscal.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dcopy : dcopy.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_dswap : dswap.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_idamax : diamax.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_idamin : diamin.o d$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_sdot : sdot.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_sasum : sasum.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_srot : srot.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_saxpy : saxpy.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_sscal : sscal.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_scopy : scopy.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_sswap : sswap.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_isamax : siamax.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_isamin : siamin.o s$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)  

test_zdot : zdot.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zasum : zasum.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_zrot : zrot.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_zaxpy : zaxpy.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zscal : zscal.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zcopy : zcopy.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_zswap : zswap.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_izamax : ziamax.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_izamin : ziamin.o z$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)

test_cdot : cdot.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_casum : casum.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_crot : crot.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)	
test_caxpy : caxpy.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_cscal : cscal.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_ccopy : ccopy.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_cswap : cswap.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_icamax : ciamax.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS)
test_icamin : ciamin.o c$(UTIL)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB) $(LIB_BLAS) $(LDFLAGS) 

clean ::
	@rm -f  *.o test_*

