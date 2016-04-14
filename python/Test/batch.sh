time CTHRES=1e-08 FLEVEL=1 FTHRES=1e-08 STHRES=1e-05 KSP="LU" mpiexec -n 2 ./solv A1.mtx b1.mtx 
time PYTHONPATH=.. mpiexec -n 2 ./solver.py --method=mumps A1.mtx b1.mtx 
