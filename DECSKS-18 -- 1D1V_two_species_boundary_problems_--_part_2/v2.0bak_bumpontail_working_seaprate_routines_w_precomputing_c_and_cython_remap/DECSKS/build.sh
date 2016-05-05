#!/bin/bash

cd lib

# build remap.so
cython -a remap.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
      -I/usr/include/python2.7 -o remap.so remap.c

cd ..
