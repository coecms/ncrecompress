ncrecompress
============

Re-compress a netcdf file in parallel using HDF5 raw IO

Uses openmp to paralleise the compression

```
export OMP_NUM_THREADS=8
pnccopy in.nc out.nc
```

The existing chunking and shuffle settings for variables are unchanged, the data is only recompressed

Todo:
 - Allow setting target compression level with argument `-d`
 - Allow setting number of threads with argument `-p`
