# -*- MakeFile -*-
# NAME: Copyright (c) 2020, Biren Patel
# LISC: MIT License
# DESC: Build unit tests for all rng modules

#------------------------------------------------------------------------------#
# Compiler Setup
#------------------------------------------------------------------------------#

cc = clang
cflag = -std=c99 -g -O3 -march=native -mavx2 -mrdrnd -m64 -pedantic -Wall \
		-Wextra -Wdouble-promotion -Wnull-dereference -Wconversion -Wcast-qual \
		-Wpacked -Wpadded

#------------------------------------------------------------------------------#
# Object Files
#------------------------------------------------------------------------------#

objects = random_test.o random_simd.o random_sisd.o random_utils.o unity.o

#------------------------------------------------------------------------------#
# Build
#------------------------------------------------------------------------------#

unit_tests.exe : $(objects)
	$(cc) $(objects) -o unit_tests.exe

unity.o : unity/unity.c unity/unity.h
	$(cc) -c unity/unity.c -o unity.o

random_test.o : random_test.c unity/unity.h timeit.h ../src/random.h
	$(cc) $(cflag) -c random_test.c -I ../src -o random_test.o

random_simd.o : ../src/random_simd.c ../src/random_simd.h ../src/random_utils.h
	$(cc) $(cflag) -c ../src/random_simd.c -o random_simd.o

random_sisd.o : ../src/random_sisd.c ../src/random_sisd.h ../src/random_utils.h \
				../src/bitarray.h
	$(cc) $(cflag) -c ../src/random_sisd.c -o random_sisd.o

random_utils.o : ../src/random_utils.c ../src/random_utils.h
	$(cc) $(cflag) -c ../src/random_utils.c -o random_utils.o

#------------------------------------------------------------------------------#
# Post-Build
#------------------------------------------------------------------------------#

clean :
	del $(objects)
