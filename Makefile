test.nc: pnccopy
	rm -f test.nc
	./pnccopy /g/data/w35/saw562/HighResMIP/EasyAerosol/1949-2015/n216e/easy_cdnc_PI-MACv2-SP_1949-2015_v2_N216e_AW_1949-2015.nc test.nc

pnccopy: pnccopy.c
	icc -g -Wall $^ -o $@ -lnetcdf -lhdf5_hl -lhdf5 -std=c99 -check-pointers=rw
