# blas-tests
Tests for checking blas correctness

## Compile

```
Simple make (default directory ../ default lib: libopenblas.a)
    make test_????
    make double
    make single
    make all

Full make: 
    make test_ddot BLAS_INC_DIR=INCLUDEDIR BLAS_LIB_DIR=LIBRARYDIR BLAS_LIBS="libopenblas.a .."

```

## Test
```
   bash call_test.sh
```