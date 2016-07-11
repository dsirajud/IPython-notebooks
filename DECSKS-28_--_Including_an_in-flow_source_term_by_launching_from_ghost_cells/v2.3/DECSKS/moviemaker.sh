#!/bin/bash

# $0 is the script name, $1 id the first ARG, $2 is second...

# cd is DECSKS/

cd plots
pwd

ffmpeg -qscale 5 -r 24 -b 9600 -i $1.png $2.mp4

mv $2.mp4 ..

# go back to DECSKS/
cd ..
