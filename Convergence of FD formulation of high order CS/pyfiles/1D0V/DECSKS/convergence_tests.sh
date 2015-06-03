#!/bin/bash

# cwd: Work/DECSKS/ , execute script from this directory!

# sys.argv[1] = (filename) simulation parameters for given run

#python main_sh.py 'params_Nx12.dat'
#echo "Nx12 complete, starting Nx24"

python main_sh.py 'params_Nx24.dat'
echo "Nx24 complete, starting Nx48"
 
python main_sh.py 'params_Nx48.dat'
echo "Nx48 complete, starting Nx96"

python main_sh.py 'params_Nx96.dat'
echo "Nx96 complete, starting Nx192"

python main_sh.py 'params_Nx192.dat'
echo "Nx192 complete, starting Nx384"

python main_sh.py 'params_Nx384.dat'
echo "Nx384 complete, next..."

python main_sh.py 'params_Nx768.dat'
echo "Nx768 complete, next..."

python main_sh.py 'params_Nx1536.dat'
echo "Nx1536 complete, done"

#python main_sh.py 'params_Nx3072.dat'
#echo "Nx3072 complete, done"
