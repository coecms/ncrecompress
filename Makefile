test: pnccopy
	py.test test.py

pnccopy: pnccopy.c
	module purge; module load netcdf/4.6.1 hdf5/1.10.2 intel-cc/2019.0.117; icc -g -Wall -Wextra $^ -o $@ -lnetcdf -lhdf5_hl -lhdf5 -std=c99 -qopenmp -check-pointers=rw -check=conversions,stack,uninit -traceback
