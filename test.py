#!/usr/bin/env python
#
# Copyright 2019 Scott Wales
#
# Author: Scott Wales <scott.wales@unimelb.edu.au>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy
import netCDF4
import subprocess
import xarray
import contextlib
import pytest
import itertools

@contextlib.contextmanager
def setup_file(path):
    ds = netCDF4.Dataset(path, mode='w')
    ds.createDimension('t',10)
    ds.createDimension('x',10)
    ds.createDimension('y',10)

    yield ds

    ds.close()

def compress_and_test(inpath, outpath):
    p = subprocess.run(['./pnccopy',inpath,outpath])
    p.check_returncode()

    a = xarray.open_dataset(inpath)
    b = xarray.open_dataset(outpath)
    
    print(a.a.encoding)
    print(b.a.encoding)

    xarray.testing.assert_identical(a,b)

def test_contiguous(tmp_path):
    inpath = tmp_path/'x.nc'
    outpath = tmp_path/'y.nc'

    with setup_file(inpath) as ds:
        v = ds.createVariable('a', float, ['t','x','y'], contiguous=True)
        v[...] = numpy.random.rand(10,10,10)

    compress_and_test(inpath, outpath)


@pytest.mark.parametrize('zlib,shuffle,fletcher32',itertools.product([True,False], repeat=3))
def test_chunk(tmp_path, zlib, shuffle, fletcher32):
    inpath = tmp_path/'x.nc'
    outpath = tmp_path/'y.nc'

    with setup_file(inpath) as ds:
        v = ds.createVariable('a', float, ['t','x','y'], chunksizes=[2,2,2], zlib=zlib, shuffle=shuffle, fletcher32=fletcher32)
        v[...] = numpy.random.rand(10,10,10)

    compress_and_test(inpath, outpath)
