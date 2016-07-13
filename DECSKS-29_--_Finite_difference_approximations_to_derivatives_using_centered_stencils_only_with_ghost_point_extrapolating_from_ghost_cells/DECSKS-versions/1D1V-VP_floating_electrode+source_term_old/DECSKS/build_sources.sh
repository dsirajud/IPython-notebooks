#!/bin/bash

cd lib

# build source_advection.so
cython -a source_advection.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o source_advection.so source_advection.c

cd ..
