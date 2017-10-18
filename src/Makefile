# This Makefile will compile all .f files in the directory
# When you add new files they will automatically be added to the list

SRC_FILES=$(wildcard *.f)

TARGET=dstup.so # Put here the name of your shared library file

$(TARGET): $(SRC_FILES)
	gfortran -shared -o $@ $^ 


david_test: david_test.f
	gfortran -o $@ $< dcopy.f dscal.f daxpy.f
