#!/bin/bash

# This script generates a directory named "KRR_binary", which is
# almost a duplicate of "KRR", except that point format is restricted
# to a dense point and the file IO is handled by LibSVM_IO_binary
# rather than LibSVM_IO. The files under "KRR_binary" should not be
# compiled with macros PointFormat and PointArrayFormat. The
# executables use the same argument list as that in the counterparts
# of "KRR". Ignore the instructions above the main function in the
# automatically generated KRR_*_binary.cpp files.

# Create directory and copy file
rm -rf KRR_binary
mkdir KRR_binary
cp KRR/* KRR_binary/
#rename .cpp _binary.cpp KRR_binary/*.cpp
#rename .hpp _binary.hpp KRR_binary/*.hpp
for file in KRR_binary/*.cpp; do
    mv $file KRR_binary/`basename $file .cpp`_binary.cpp
done
for file in KRR_binary/*.hpp; do
    mv $file KRR_binary/`basename $file .hpp`_binary.hpp
done

# Modification in KRR_Common_binary.hpp
file="KRR_binary/KRR_Common_binary.hpp"
sed -i'' -e 's/KRR_COMMON/KRR_COMMON_BINARY/g' $file

# Modification in KRR_Common_binary.cpp
file="KRR_binary/KRR_Common_binary.cpp"
sed -i'' -e 's/KRR_Common.hpp/KRR_Common_binary.hpp/g' $file

# Modification in other KRR_*_binary.cpp
method_list=("Standard" "Nystrom" "Fourier" "RLCM" "BlockDiag")
for method in "${method_list[@]}"; do
    file="KRR_binary/KRR_${method}_binary.cpp"
    sed -i'' -e 's/KRR_Common.hpp/KRR_Common_binary.hpp/g' $file
    sed -i'' -e 's/PointArrayFormat/DPointArray/g' $file
    sed -i'' -e 's/PointFormat/DPoint/g' $file
    sed -i'' -e 's/LibSVM_IO/LibSVM_IO_binary/g' $file
done

# Modification in Makefile
file="KRR_binary/Makefile"
sed -i'' -e 's/KRR_Common.o/KRR_Common_binary.o/g' $file
sed -i'' -e 's/.cpp/_binary.cpp/g' $file
sed -i'' -e 's/%_binary.cpp/%.cpp/g' $file
sed -i'' -e 's/-D PointFormat=DPoint -D PointArrayFormat=DPointArray//g' $file
sed -i'' -e 's/_DPoint/_binary/g' $file
sed -i'' -e '/SPoint/d' $file
