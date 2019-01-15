#!/bin/bash

# The user must edit this part ------------------------------------------------

# Identify the operating system. Either LINUX or MAC.
export SYS=MAC

# Set C++ compiler. Examples: g++ (Linux), g++-8 (Mac).
export CXX=g++-8

# Set the use of OPENMP. The C++ compiler typically comes with OPENMP
# support.
export USE_OPENMP=1

# The LAPACK support, if any, is often single threaded. One may use
# one (and only one) of the following multithreaded libraries for
# substitution: OPENBLAS, ESSL, MKL, and ACCELERATE.

# Set the use of OPENBLAS. If USE_OPENBLAS is zero, the other related
# variables are ignored.
export USE_OPENBLAS=0
export OPENBLAS_LIBS="-lopenblaso"
export OPENBLAS_INCL_PATH=
export OPENBLAS_LIB_PATH=

# Set the use of ESSL. If USE_ESSL is zero, the other related
# variables are ignored.
export USE_ESSL=0
export ESSL_CXXFLAGS="-q64"
export ESSL_LIBS="-lesslsmp -lxlf90_r -lxl -lxlsmp -lxlfmath -llapack"
export ESSL_INCL_PATH="-I/usr/include/"
export ESSL_LIB_PATH="-L/usr/lib/ -L/opt/ibm/xlf/15.1.4/lib/ -L/opt/ibm/xlsmp/4.1.2/lib/"

# Set the use of MKL. If USE_MKL is zero, the other related variables
# are ignored.
export USE_MKL=0
export MKL_CXXFLAGS="-openmp -mkl"
export MKL_LIBS="-openmp -mkl"

# Set the use of ACCELERATE.
export USE_ACCELERATE=1

# Set the use of 64-bit integers. ALWAYS SET IT TO ZERO.
#
# The current support for USE_LONG is bogus and is not fully
# functional. The situation is delicate and the root cause lies in the
# 64-bit integer support of LAPACK. Long story short, not all LAPACK
# libraries fully support 64-bit.
export USE_LONG=0

# Set the use of METIS. ALWAYS SET IT TO ZERO.
#
# The current version of the library does not require METIS.
export HAS_METIS=0
export METIS_INCL_PATH=
export METIS_LIB_PATH=

# Set the use of Matern kernel.
export USE_MATERN=1

# If the Matern kernel is used, one must also provide the Fortran
# compiler. Examples: gfortran (Linux), gfortran-8 (Mac).
export FC=gfortran-8

# The user stops editing here -------------------------------------------------

export CMATRIX_DIR=$PWD

export OUTFILE=$CMATRIX_DIR/Makefile
cat $OUTFILE.in                                                 >  $OUTFILE

export OUTFILE=$CMATRIX_DIR/src/Makefile
cat $OUTFILE.in                                                 >  $OUTFILE

export SRC_SUBDIRS=( Misc Solvers Kernels TestFunctions Matrices KRR GP )

for DIR in "${SRC_SUBDIRS[@]}"; do
    export OUTFILE=$CMATRIX_DIR/src/$DIR/Makefile
    echo "SYS = $SYS"                                           >  $OUTFILE
    echo "CXX = $CXX"                                           >> $OUTFILE
    echo "USE_OPENMP = $USE_OPENMP"                             >> $OUTFILE
    echo "USE_OPENBLAS = $USE_OPENBLAS"                         >> $OUTFILE
    echo "OPENBLAS_INCL_PATH = $OPENBLAS_INCL_PATH"             >> $OUTFILE
    echo "USE_ESSL = $USE_ESSL"                                 >> $OUTFILE
    echo "ESSL_CXXFLAGS = $ESSL_CXXFLAGS"                       >> $OUTFILE
    echo "ESSL_INCL_PATH = $ESSL_INCL_PATH"                     >> $OUTFILE
    echo "USE_MKL = $USE_MKL"                                   >> $OUTFILE
    echo "MKL_CXXFLAGS = $MKL_CXXFLAGS"                         >> $OUTFILE
    echo "USE_ACCELERATE = $USE_ACCELERATE"                     >> $OUTFILE
    echo "USE_LONG = $USE_LONG"                                 >> $OUTFILE
    echo "HAS_METIS = $HAS_METIS"                               >> $OUTFILE
    echo "METIS_INCL_PATH = $METIS_INCL_PATH"                   >> $OUTFILE
    echo "USE_MATERN = $USE_MATERN"                             >> $OUTFILE
    echo "FC = $FC"                                             >> $OUTFILE
    cat $OUTFILE.in                                             >> $OUTFILE
done

export OUTFILE=$CMATRIX_DIR/test/Makefile
echo "CMATRIX_DIR = $CMATRIX_DIR"                               >  $OUTFILE
echo "SYS = $SYS"                                               >> $OUTFILE
echo "CXX = $CXX"                                               >> $OUTFILE
echo "USE_OPENMP = $USE_OPENMP"                                 >> $OUTFILE
echo "USE_OPENBLAS = $USE_OPENBLAS"                             >> $OUTFILE
echo "OPENBLAS_LIBS = $OPENBLAS_LIBS"                           >> $OUTFILE
echo "OPENBLAS_INCL_PATH = $OPENBLAS_INCL_PATH"                 >> $OUTFILE
echo "OPENBLAS_LIB_PATH = $OPENBLAS_LIB_PATH"                   >> $OUTFILE
echo "USE_ESSL = $USE_ESSL"                                     >> $OUTFILE
echo "ESSL_CXXFLAGS = $ESSL_CXXFLAGS"                           >> $OUTFILE
echo "ESSL_LIBS = $ESSL_LIBS"                                   >> $OUTFILE
echo "ESSL_INCL_PATH = $ESSL_INCL_PATH"                         >> $OUTFILE
echo "ESSL_LIB_PATH = $ESSL_LIB_PATH"                           >> $OUTFILE
echo "USE_MKL = $USE_MKL"                                       >> $OUTFILE
echo "MKL_CXXFLAGS = $MKL_CXXFLAGS"                             >> $OUTFILE
echo "MKL_LIBS = $MKL_LIBS"                                     >> $OUTFILE
echo "USE_ACCELERATE = $USE_ACCELERATE"                         >> $OUTFILE
echo "USE_LONG = $USE_LONG"                                     >> $OUTFILE
echo "HAS_METIS = $HAS_METIS"                                   >> $OUTFILE
echo "METIS_INCL_PATH = $METIS_INCL_PATH"                       >> $OUTFILE
echo "METIS_LIB_PATH = $METIS_LIB_PATH"                         >> $OUTFILE
echo "USE_MATERN = $USE_MATERN"                                 >> $OUTFILE
echo "FC = $FC"                                                 >> $OUTFILE
cat $OUTFILE.in                                                 >> $OUTFILE

export OUTFILE=$CMATRIX_DIR/app/Makefile
cat $OUTFILE.in                                                 >  $OUTFILE

export APP_SUBDIRS=( KRR GP )

for DIR in "${APP_SUBDIRS[@]}"; do
    export OUTFILE=$CMATRIX_DIR/app/$DIR/Makefile
    echo "CMATRIX_DIR = $CMATRIX_DIR"                           >  $OUTFILE
    echo "SYS = $SYS"                                           >> $OUTFILE
    echo "CXX = $CXX"                                           >> $OUTFILE
    echo "USE_OPENMP = $USE_OPENMP"                             >> $OUTFILE
    echo "USE_OPENBLAS = $USE_OPENBLAS"                         >> $OUTFILE
    echo "OPENBLAS_LIBS = $OPENBLAS_LIBS"                       >> $OUTFILE
    echo "OPENBLAS_INCL_PATH = $OPENBLAS_INCL_PATH"             >> $OUTFILE
    echo "OPENBLAS_LIB_PATH = $OPENBLAS_LIB_PATH"               >> $OUTFILE
    echo "USE_ESSL = $USE_ESSL"                                 >> $OUTFILE
    echo "ESSL_CXXFLAGS = $ESSL_CXXFLAGS"                       >> $OUTFILE
    echo "ESSL_LIBS = $ESSL_LIBS"                               >> $OUTFILE
    echo "ESSL_INCL_PATH = $ESSL_INCL_PATH"                     >> $OUTFILE
    echo "ESSL_LIB_PATH = $ESSL_LIB_PATH"                       >> $OUTFILE
    echo "USE_MKL = $USE_MKL"                                   >> $OUTFILE
    echo "MKL_CXXFLAGS = $MKL_CXXFLAGS"                         >> $OUTFILE
    echo "MKL_LIBS = $MKL_LIBS"                                 >> $OUTFILE
    echo "USE_ACCELERATE = $USE_ACCELERATE"                     >> $OUTFILE
    echo "USE_LONG = $USE_LONG"                                 >> $OUTFILE
    echo "HAS_METIS = $HAS_METIS"                               >> $OUTFILE
    echo "METIS_INCL_PATH = $METIS_INCL_PATH"                   >> $OUTFILE
    echo "METIS_LIB_PATH = $METIS_LIB_PATH"                     >> $OUTFILE
    echo "USE_MATERN = $USE_MATERN"                             >> $OUTFILE
    echo "FC = $FC"                                             >> $OUTFILE
    cat $OUTFILE.in                                             >> $OUTFILE
done
