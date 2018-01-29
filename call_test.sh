#!/bin/bash
for x in $(ls test*)
do
echo "Tests $x"
export OPENBLAS_INCX=1
export OPENBLAS_INCY=1
echo "INCX=${OPENBLAS_INCX} INCY=${OPENBLAS_INCY}"
./$x 5  177 
export OPENBLAS_INCY=4
echo "INCX=${OPENBLAS_INCX} INCY=${OPENBLAS_INCY}"

./$x 5 177
export OPENBLAS_INCX=4
echo "INCX=${OPENBLAS_INCX} INCY=${OPENBLAS_INCY}"

./$x 5 177
export OPENBLAS_INCX=7
export OPENBLAS_INCY=7
echo "INCX=${OPENBLAS_INCX} INCY=${OPENBLAS_INCY}"

./$x 5 177
done
