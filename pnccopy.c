/** 
 * Copyright 2019 Scott Wales
 *
 * \author  Scott Wales <scott.wales@unimelb.edu.au>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "netcdf.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#define ERR_IF(x) err_if_(x, #x, __FILE__, __LINE__)
void err_if_(bool err, const char * message, const char * file, size_t line) {
    if (err) {
        fprintf(stderr, "%s:%zu %s\n", file, line, message);
        exit(-1);
    }
}

#define NC_ERR(x) nc_err_(x, #x, __FILE__, __LINE__)
void nc_err_(int err, const char * message, const char * file, size_t line) {
    const char * nc_str = nc_strerror(err);
    err_if_(err != NC_NOERR, nc_str, file, line);
}

#define H5_ERR(x) h5_err_(x, #x, __FILE__, __LINE__)
hid_t h5_err_(hid_t err, const char * message, const char * file, size_t line) {
    if (err < 0) {
        H5Eprint(H5Eget_current_stack(), stderr);
        err_if_(err < 0, message, file, line);
    }
    return err;
}


void copy_dim(int dimid, int ncidin, int ncidout) {
    char name[NC_MAX_NAME+1];
    size_t len;
    int dimidout;

    NC_ERR(nc_inq_dim(ncidin, dimid, name, &len));
    NC_ERR(nc_def_dim(ncidout, name, len, &dimidout));

    ERR_IF(dimid != dimidout);
}

void copy_att(int varid, int attid, int ncidin, int ncidout) {
    char name[NC_MAX_NAME+1];
    nc_type xtype;
    size_t len;

    NC_ERR(nc_inq_attname(ncidin, varid, attid, name));
    NC_ERR(nc_inq_att(ncidin, varid, name, &xtype, &len));

    char buffer[len * 8];
    NC_ERR(nc_get_att(ncidin, varid, name, buffer));

    NC_ERR(nc_put_att(ncidout, varid, name, xtype, len, buffer));
}

void copy_var(int varid, int ncidin, int ncidout) {
    char name[NC_MAX_NAME+1];
    int ndims;
    int natts;
    nc_type xtype;
    int varidout;

    NC_ERR(nc_inq_var(ncidin, varid, name, &xtype, &ndims, NULL, &natts));

    int dimids[ndims];
    int shuffle;
    int deflate;
    int deflate_level;
    int contiguous;
    size_t chunksizes[ndims];

    NC_ERR(nc_inq_var(ncidin, varid, NULL, NULL, NULL, dimids, NULL));
    NC_ERR(nc_inq_var_deflate(ncidin, varid, &shuffle, &deflate, &deflate_level));
    NC_ERR(nc_inq_var_chunking(ncidin, varid, &contiguous, chunksizes));

    NC_ERR(nc_def_var(ncidout, name, xtype, ndims, dimids, &varidout));
    NC_ERR(nc_def_var_deflate(ncidout, varidout, shuffle, deflate, deflate_level));
    NC_ERR(nc_def_var_chunking(ncidout, varidout, contiguous, chunksizes));

    ERR_IF(varid != varidout);

    for (int a=0; a<natts; ++a) {
        copy_att(varid, a, ncidin, ncidout);
    }
}

void copy_structure(const char * pathin, const char * pathout) {
    
    int ncidin;
    NC_ERR(nc_open(pathin, NC_NOWRITE, &ncidin));

    int ncidout;
    NC_ERR(nc_create(pathout, NC_NOCLOBBER | NC_NETCDF4 | NC_CLASSIC_MODEL, &ncidout));

    int ndims, nvars, natts;
    NC_ERR(nc_inq(ncidin, &ndims, &nvars, &natts, NULL));

    printf("Dims %d, Vars %d, Atts %d\n", ndims, nvars, natts);

    for (int d=0; d<ndims; ++d) {
        copy_dim(d, ncidin, ncidout);
    }

    for (int v=0; v<nvars; ++v) {
        copy_var(v, ncidin, ncidout);
    }

    for (int a=0; a<natts; ++a) {
        copy_att(NC_GLOBAL, a, ncidin, ncidout);
    }

    NC_ERR(nc_close(ncidout));
    NC_ERR(nc_close(ncidin));
}

void copy_chunks(int ndims, const hsize_t size[], const hsize_t chunksize[], hid_t data_in, hid_t data_out) {
    size_t nchunks = 1;
    size_t nchunks_d[ndims];
    for (int d=0; d<ndims; ++d) {
        nchunks_d[d] = size[d] / chunksize[d];
        nchunks *= nchunks_d[d];
    }

    hid_t pread = H5_ERR(H5Pcreate(H5P_DATATYPE_ACCESS));

    hsize_t buffer_size = 0;
    char * buffer = NULL;

    for (size_t c=0; c<nchunks; ++c) {
        hsize_t offset[ndims];
        size_t tmp = c;
        for (int d=0; d<ndims; ++d) {
            offset[d] = tmp % nchunks_d[d] * chunksize[d];
            tmp = tmp / nchunks_d[d];
        }

        hsize_t chunk_size;
        H5_ERR(H5Dget_chunk_storage_size(data_in, offset, &chunk_size));

        if (chunk_size < buffer_size) {
            buffer_size = chunk_size;
            buffer = realloc(buffer, buffer_size);
        }

        uint32_t filter_mask = 0;
        H5_ERR(H5DOread_chunk(data_in, H5P_DEFAULT, offset, &filter_mask, buffer));
    }

    H5_ERR(H5Pclose(pread));
}

herr_t copy_object(hid_t oid, const char * name, const H5O_info_t * info, void * op_data) {
    if (info->type != H5O_TYPE_DATASET) {return 0;}

    printf("%s %d\n", name, info->type);

    hid_t data_in = H5_ERR(H5Dopen(oid, name, H5P_DEFAULT));


    // Get chunking info
    hid_t pcreate = H5_ERR(H5Dget_create_plist(data_in));
    H5D_layout_t layout = H5_ERR(H5Pget_layout(pcreate));
    
    if (layout == H5D_CHUNKED) {
        hid_t space = H5_ERR(H5Dget_space(data_in));

        int ndims = H5_ERR(H5Sget_simple_extent_ndims(space));
        hsize_t size[ndims];
        H5_ERR(H5Sget_simple_extent_dims(space, size, NULL));

        hsize_t chunksize[ndims];
        H5_ERR(H5Pget_chunk(pcreate, ndims, chunksize));

        hid_t file_out = *(hid_t*)op_data;
        hid_t data_out = H5_ERR(H5Dopen(file_out, name, H5P_DEFAULT));

        copy_chunks(ndims, size, chunksize, data_in, data_out);

        H5_ERR(H5Sclose(space));
        H5_ERR(H5Dclose(data_out));
    }

    H5_ERR(H5Pclose(pcreate));
    H5_ERR(H5Dclose(data_in));
    return 0;
}

void copy_data(const char * pathin, const char * pathout) {
    hid_t file_in = H5_ERR(H5Fopen(pathin, H5F_ACC_RDONLY, H5P_DEFAULT));
    hid_t file_out = H5_ERR(H5Fopen(pathout, H5F_ACC_RDWR, H5P_DEFAULT));

    H5_ERR(H5Ovisit(file_in, H5_INDEX_NAME, H5_ITER_NATIVE, copy_object, &file_out));

    H5_ERR(H5Fclose(file_out));
    H5_ERR(H5Fclose(file_in));
}

int main(int argc, char ** argv) {
    copy_structure(argv[1], argv[2]);
    copy_data(argv[1], argv[2]);
}
