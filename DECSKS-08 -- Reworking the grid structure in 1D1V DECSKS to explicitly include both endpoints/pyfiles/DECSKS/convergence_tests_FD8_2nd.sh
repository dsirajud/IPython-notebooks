#!/bin/bash

# cwd: Work/DECSKS/ , execute script from this directory!

# sys.argv[1] = (filename) simulation parameters for given run

python main_sh.py 'params_Nx9.dat'
echo "Nx9 complete, starting Nx18"
 
python main_sh.py 'params_Nx18.dat'
echo "Nx18 complete, starting Nx36"

python main_sh.py 'params_Nx36.dat'
echo "Nx36 complete, starting Nx72"

python main_sh.py 'params_Nx72.dat'
echo "Nx72 complete, starting Nx144"

python main_sh.py 'params_Nx144.dat'
echo "Nx144 complete, next..."

python main_sh.py 'params_Nx288.dat'
echo "Nx288 complete, next..."

python main_sh.py 'params_Nx576.dat'
echo "Nx576 complete, done"

