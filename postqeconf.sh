#!/bin/bash 

./configure  DFLAGS="-D__FFTW" CFLAGS=-fPIC FFLAGS="-g -fPIC" --disable-parallel
