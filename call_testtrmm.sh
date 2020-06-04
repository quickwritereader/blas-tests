#!/bin/bash

k_start=1
k_end=480
k_blocksize=1

# armie -msve-vector-bits=256 -i libinscount_emulated.so -- ./test_bl_dtrmm.x     $k $k
# ~/qemu.git/aarch64-linux-user/qemu-aarch64 -cpu max,sve256=on ./test_bl_dtrmm.x  $k $k
CMD_PRFX="" #prefix 
LOGERR_FILE="logerr.t"
LOGOUT_FILE="logout.t"

rm -rf $LOGERR_FILE $LOGOUT_FILE

for (( k=k_start; k<=k_end; k+=k_blocksize ))
do

BLAS_SIDE='L' BLAS_UPLO='U' BLAS_TRANS='N' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='L' BLAS_UPLO='U' BLAS_TRANS='N' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE

BLAS_SIDE='L' BLAS_UPLO='U' BLAS_TRANS='T' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='L' BLAS_UPLO='U' BLAS_TRANS='T' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE

BLAS_SIDE='L' BLAS_UPLO='L' BLAS_TRANS='N' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='L' BLAS_UPLO='L' BLAS_TRANS='N' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE

BLAS_SIDE='L' BLAS_UPLO='L' BLAS_TRANS='T' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='L' BLAS_UPLO='L' BLAS_TRANS='T' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE

BLAS_SIDE='R' BLAS_UPLO='U' BLAS_TRANS='N' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='R' BLAS_UPLO='U' BLAS_TRANS='N' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE

BLAS_SIDE='R' BLAS_UPLO='U' BLAS_TRANS='T' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='R' BLAS_UPLO='U' BLAS_TRANS='T' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE

BLAS_SIDE='R' BLAS_UPLO='L' BLAS_TRANS='N' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='R' BLAS_UPLO='L' BLAS_TRANS='N' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE

BLAS_SIDE='R' BLAS_UPLO='L' BLAS_TRANS='T' BLAS_DIAG='U' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE
BLAS_SIDE='R' BLAS_UPLO='L' BLAS_TRANS='T' BLAS_DIAG='N' $CMD_PRFX ./test_dtrmm $k $k 2>>$LOGERR_FILE 1>>$LOGOUT_FILE


done

echo "END" >>$LOGERR_FILE
echo "END"
