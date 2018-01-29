#!/bin/bash
for x in $(ls test*)
do
echo "Tests $x"
export OPENBLAS_INCX=1
export OPENBLAS_INCY=1
./$x 8  377 
export OPENBLAS_INCY=4
./$x 8 377
export OPENBLAS_INCX=4
./$x 8 377
export OPENBLAS_INCX=7
export OPENBLAS_INCY=7
./$x 8 377
done
