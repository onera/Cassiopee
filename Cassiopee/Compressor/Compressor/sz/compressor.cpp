/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <array>
#include <iostream>
#define K_ARRAY_UNIQUE_SYMBOL
#include "Python.h"
#include "kcore.h"
#include "sz_omp.h"
#include "zlib.h"

namespace sz {
enum sz_ikeys {
    snapShotCmprStep = 0,
    withLinearRegression,
    protectValueRange,
    sampleDistance,
    quantization_intervals,
    max_quant_intervals,
    predThreshold,
    szMode,
    losslessCompressor,
    gzipMode,
    zstdMode,
    randomAccess,
    errorBoundMode,
    absErrBound,
    relBoundRatio,
    psnr,
    pw_relBoundRatio,
    accelerate_pw_rel_compression,
    end_guard
};
/** Définition des clefs pour SZ
 */
constexpr std::array<const char *, end_guard> sz_keys = {"snapshotCmprStep",
                                                         "withLinearRegression",
                                                         "protectValueRange",
                                                         "sampleDistance",
                                                         "quantization_intervals",
                                                         "max_quant_intervals",
                                                         "predThreshold",
                                                         "szMode",
                                                         "losslessCompressor",
                                                         "gzipMode",
                                                         "zstdMode",
                                                         "randomAccess",
                                                         "errorBoundMode",
                                                         "absErrBound",
                                                         "relBoundRatio",
                                                         "psnr",
                                                         "pw_relBoundRatio",
                                                         "accelerate_pw_rel_compression"};

void
init_parameters(sz_params &parameters)
{
    /*sz_params sz;
    memset(&sz, 0, sizeof(sz_params));
    sz.sol_ID = SZ;
    sz.sampleDistance = 100;
    sz.quantization_intervals = 0;
    sz.max_quant_intervals = 65536;
    sz.predThreshold = 0.98;
    sz.szMode = SZ_BEST_COMPRESSION;
    sz.losslessCompressor = ZSTD_COMPRESSOR;
    sz.gzipMode = 1;
    sz.errorBoundMode = REL;
    sz.absErrBound = 1E-6;
    sz.relBoundRatio = 1E-5;

    SZ_Init_Params(&sz);*/
    memset(&parameters, 0, sizeof(sz_params));
    parameters.sol_ID = SZ;
    parameters.protectValueRange      = 1;
    parameters.sampleDistance         = 100;
    parameters.quantization_intervals = 0;
    parameters.max_quant_intervals    = 65536;
    parameters.predThreshold          = 0.98;
    parameters.szMode                 = SZ_BEST_COMPRESSION;
    parameters.losslessCompressor     = ZSTD_COMPRESSOR;
    parameters.gzipMode               = 1;//Fastest compression or 19;// <= Best compression for zstd
    parameters.errorBoundMode         = REL;
    parameters.absErrBound            = 1.E-6;
    parameters.relBoundRatio          = 1.E-5;
    parameters.losslessCompressor     = ZSTD_COMPRESSOR;
}

void display_parameters(const sz_params &params)
{
    if(params.sol_ID == SZ)
        printf("compressor Name:        \t SZ\n");
    else if(params.sol_ID == SZ_Transpose)
        printf("compressor Name:        \t SZ_Transpose\n");
    else
        printf("compressor Name:        \t Other compressor\n");
    switch(params.dataType)
    {
    case SZ_FLOAT:
        printf("Data type:                      \t FLOAT\n");
        printf("min value of raw data:          \t %f\n", params.fmin);
        printf("max value of raw data:          \t %f\n", params.fmax);        
        break;
    case SZ_DOUBLE:
        printf("Data type:                      \t DOUBLE\n");
        printf("min value of raw data:          \t %f\n", params.dmin);
        printf("max value of raw data:          \t %f\n", params.dmax);    
        break;
    case SZ_INT8:
        printf("Data type:                      \t INT8\n");
        break;  
    case SZ_INT16:
        printf("Data type:                      \t INT16\n");
        break;
    case SZ_INT32:
        printf("Data type:                      \t INT32\n");
        break;  
    case SZ_INT64:
        printf("Data type:                      \t INT64\n");
        break;  
    case SZ_UINT8:
        printf("Data type:                      \t UINT8\n");
        break;  
    case SZ_UINT16:
        printf("Data type:                      \t UINT16\n");
        break;
    case SZ_UINT32:
        printf("Data type:                      \t UINT32\n");
        break;  
    case SZ_UINT64:
        printf("Data type:                      \t UINT64\n");
        break;              
    }
        
    printf("dataEndianType (prior raw data):\t %s\n", dataEndianType==BIG_ENDIAN_DATA?"BIG_ENDIAN":"LITTLE_ENDIAN");
    printf("sysEndianType (at compression): \t %s\n", sysEndianType==1?"BIG_ENDIAN":"LITTLE_ENDIAN");
    printf("sampleDistance:                 \t %d\n", params.sampleDistance);
    printf("predThreshold:                  \t %f\n", params.predThreshold);
    switch(params.szMode)
    {
    case SZ_BEST_SPEED:
        printf("szMode:                         \t SZ_BEST_SPEED (without Gzip)\n");
        break;
    case SZ_BEST_COMPRESSION:
        printf("szMode:                         \t SZ_BEST_COMPRESSION (with Zstd or Gzip)\n");
        break;
    }
    switch(params.gzipMode)
    {
    case Z_BEST_SPEED:
        printf("gzipMode:                       \t Z_BEST_SPEED\n");
        break;
    case Z_DEFAULT_COMPRESSION:
        printf("gzipMode:                       \t Z_BEST_SPEED\n");
        break;  
    case Z_BEST_COMPRESSION:
        printf("gzipMode:                       \t Z_BEST_COMPRESSION\n");
        break;
    }
    
    switch(params.errorBoundMode)
    {
    case ABS:
        printf("errBoundMode:                   \t ABS\n");
        printf("absErrBound:                    \t %g\n", params.absErrBound);
        break;
    case REL:
        printf("errBoundMode:                   \t REL (based on value_range extent)\n");
        printf("relBoundRatio:                  \t %g\n", params.relBoundRatio);
        break;
    case ABS_AND_REL:
        printf("errBoundMode:                   \t ABS_AND_REL\n");
        printf("absErrBound:                    \t %g\n", params.absErrBound);
        printf("relBoundRatio:                  \t %g\n", params.relBoundRatio);
        break;
    case ABS_OR_REL:
        printf("errBoundMode:                   \t ABS_OR_REL\n");
        printf("absErrBound:                    \t %g\n", params.absErrBound);
        printf("relBoundRatio:                  \t %g\n", params.relBoundRatio);
        break;
    case PSNR:
        printf("errBoundMode:                   \t PSNR\n");
        printf("psnr:                           \t %g\n", params.psnr);
        break;
    case PW_REL:
        printf("errBoundMode:                   \t PW_REL\n");
        break;
    case ABS_AND_PW_REL:
        printf("errBoundMode:                   \t ABS_AND_PW_REL\n");
        printf("absErrBound:                    \t %f\n", params.absErrBound);
        break;
    case ABS_OR_PW_REL:
        printf("errBoundMode:                   \t ABS_OR_PW_REL\n");
        printf("absErrBound:                    \t %f\n", params.absErrBound);
        break;
    case REL_AND_PW_REL:
        printf("errBoundMode:                   \t REL_AND_PW_REL\n");
        printf("range_relBoundRatio:            \t %f\n", params.relBoundRatio);
        break;
    case REL_OR_PW_REL:
        printf("errBoundMode:                   \t REL_OR_PW_REL\n");
        printf("range_relBoundRatio:            \t %f\n", params.relBoundRatio);
        break;
    }
    
    if(params.errorBoundMode>=PW_REL && params.errorBoundMode<=REL_OR_PW_REL)
    {
        printf("pw_relBoundRatio:               \t %f\n", params.pw_relBoundRatio);
        //printf("segment_size:                   \t %d\n", params.segment_size);
        switch(params.pwr_type)
        {
        case SZ_PWR_MIN_TYPE:
            printf("pwrType:                    \t SZ_PWR_MIN_TYPE\n");
            break;
        case SZ_PWR_AVG_TYPE:
            printf("pwrType:                    \t SZ_PWR_AVG_TYPE\n");
            break;
        case SZ_PWR_MAX_TYPE:
            printf("pwrType:                    \t SZ_PWR_MAX_TYPE\n");
            break;
        }
    }
    fflush(stdout);
}

bool
set_parameters_from_dictionnary(PyObject *options, sz_params &parameters)
{
    PyObject * key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(options, &pos, &key, &value)) {
#if PY_VERSION_HEX >= 0x03000000
        std::string ckey = PyUnicode_AsUTF8(key);
#else
        std::string ckey = PyString_AsString(key);
#endif
        if (ckey == sz_keys[quantization_intervals]) {
            if (!PyLong_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an integer type for quantization intervals");
                return false;
            }
            long val                          = PyLong_AsLong(value);
            parameters.quantization_intervals = val;
        }
        else if (ckey == sz_keys[max_quant_intervals]) {
            if (!PyLong_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an integer type for max quant intervals");
                return false;
            }
            long val                       = PyLong_AsLong(value);
            parameters.max_quant_intervals = val;
        }
        else if (ckey == sz_keys[sampleDistance]) {
            if (!PyLong_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an integer type for sample distance");
                return false;
            }
            long val                  = PyLong_AsLong(value);
            parameters.sampleDistance = val;
        }
        else if (ckey == sz_keys[predThreshold]) {
            if (!PyFloat_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an float type for pred threshold");
                return false;
            }
            double val               = PyFloat_AsDouble(value);
            parameters.predThreshold = val;
        }
        else if (ckey == sz_keys[szMode]) {
            if (!PyLong_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an integer type for sz mode");
                return false;
            }
            long val          = PyLong_AsLong(value);
            parameters.szMode = val;
        }
        else if (ckey == sz_keys[gzipMode]) {
            if (!PyLong_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an integer type for gzip mode");
                return false;
            }
            long val            = PyLong_AsLong(value);
            parameters.gzipMode = val;
        }
        else if (ckey == sz_keys[errorBoundMode]) {
            if (!PyLong_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an integer type for error bound mode");
                return false;
            }
            long val                  = PyLong_AsLong(value);
            parameters.errorBoundMode = val;
        }
        else if (ckey == sz_keys[absErrBound]) {
            if (!PyFloat_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an float type for abs error bound");
                return false;
            }
            double val             = PyFloat_AsDouble(value);
            parameters.absErrBound = val;
        }
        else if (ckey == sz_keys[relBoundRatio]) {
            if (!PyFloat_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an float type for rel bound ratio");
                return false;
            }
            double val               = PyFloat_AsDouble(value);
            parameters.relBoundRatio = val;
        }
        else if (ckey == sz_keys[pw_relBoundRatio]) {
            if (!PyFloat_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an float type for point wise rel bound ratio");
                return false;
            }
            double val                  = PyFloat_AsDouble(value);
            parameters.pw_relBoundRatio = val;
        }
        else if (ckey == sz_keys[psnr]) {
            if (!PyLong_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting an integer type for psnr");
                return false;
            }
            long val        = PyLong_AsLong(value);
            parameters.psnr = val;
        }
        else if (ckey == sz_keys[protectValueRange]) {
            if (!PyBool_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting a boolean type for protectValueRange");
                return false;
            }
            parameters.protectValueRange = (PyObject_IsTrue(value) ? 1 : 0);
        }
        else if (ckey == sz_keys[accelerate_pw_rel_compression]) {
            if (!PyBool_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting a boolean type for accelerate_pw_rel_compression");
                return false;
            }
            parameters.accelerate_pw_rel_compression = (PyObject_IsTrue(value) ? 1 : 0);
        }
        else if (ckey == sz_keys[randomAccess]) {
            if (!PyBool_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting a boolean type for protectValueRange");
                return false;
            }
            parameters.randomAccess = (PyObject_IsTrue(value) ? 1 : 0);
        }
        else if (ckey == sz_keys[withLinearRegression]) {
            if (!PyBool_Check(value)) {
                PyErr_SetString(PyExc_TypeError, "Waiting a boolean type for withLinearFegression");
                return false;
            }
            parameters.withRegression = (PyObject_IsTrue(value) ? SZ_WITH_LINEAR_REGRESSION : SZ_NO_REGRESSION );
        }
        else {
            PyErr_SetString(PyExc_Warning, "Key for sz options doesn't exist !");
            return false;
        }
    }
    return true;
}

const char *compress_doc = R"RAW( 
        Compress an array or a list of arrays of double values. 
        ------------------------------------------------------
        Usage : cpr = sz.pack(array[, options])
              or
                lst_cpr = sz.pack([arr_1,arr_2,...,arr_n][, options])
        with options could be :
            - a file describing the options to compress the double dataset
            - a dictionnary with options to change

        Some options are possible to have better compression (sorted by category) :
            [System environment options <= no need to change unless necessary]
             ----------------------------------------------------------------
            - endianType : either LITTLE_ENDIAN_DATA or BIG_ENDIAN_DATA
                           x86, x64 and arm adopt LITTLE_ENDIAN_DATA
                           PowerPC (PPC), MAC OS, and KEIL C51 adopt BIG_ENDIAN_DATA

            [Compression parameters <= no need to change unless necessary]
             ------------------------------------------------------------
            - snapshotCmprStep : is used to define the period of spatial-compression during the time-based compression 
              (In order to support time-based compression, you need to enable time-based compression by using --enable-timecmpr during the compilation.)

            - withLinearRegression : Yes means using SZ 1.4, No means using SZ 2.1 (default YES)

            - protectValueRange : allows to preserve the value range for the decompressed data (value: YES or NO, default YES)
                                  Switching on this option may introduce a little execution cost in decompression, but no impact to compression time at all.

            - sampleDistance : determins the number of samples used to optimize the # quantization intervals
                               For example, sampleDistance=50 means 1/50=2% of data points are sample points.
                               Default value is 100

            - quantization_intervals : The number of quantization intervals should be always set to an "even" number!
                                       If it is set to 0, SZ will autamatically search for an optimized setting.
                   :  it has be to no less than 4 and no greater than 65536, such as 256.
                                       Default value is 0.

            - max_quant_intervals : maximum quantization interval is valid only when quantization_intervals=0 (i.e., let the sz compressor optimize the intervals)
                                    In general, this setting does not change the compression ratio/rate, but only affect the compression speed to a certain extent (only 10% in general).
                                    The high values of max_quant_intervals, the lower compression speed, but the higher ability the compressor can reach high compression ratios for high-precision compression.
                                    As for low-precision compression (i.e., high error bound such as 1E-2), max_quant_intervals could be set to 256 or 65536.                                    #As for pretty-high-precision demand (i.e., fairly small error bound such as 1E-6), max_quant_intervals could be set to 2097152(=2^21).
                                    Default value is 65536

            - predThreshold: the threshold to determine the ratio of predictable data over all data
                             predThreshold = 0.97 means 97% of data will be predictable.
                             The default value is 0.99

            - SZ_Mode : two options: SZ_BEST_SPEED or SZ_BEST_COMPRESSION (default value is SZ_BEST_COMPRESSION)

            - losslessCompressor : Select the lossless compression techniques after the lossy compression: either ZSTD_COMPRESSOR or GZIP_COMPRSSOR
                                   Default value is ZSTD_COMPRESSOR

            - GZIP_Mode (only valid when losslessCompressor is GZIP_COMPRESSOR) : 
                         Gzip_NO_COMPRESSION, Gzip_DEFAULT_COMPRESSION, Gzip_BEST_SPEED or Gzip_BEST_COMPRESSION (default Gzip_BEST_SPEED)
                         Note: this parameter setting is valid only if szMode = SZ_BEST_COMPRESION.

            - zstdMode (only valid when losslessCompressor is GZIP_COMPRESSOR or ZSTD_COMPRESSOR) :
                       If losslessCompressor = ZSTD_COMPRESSOR, there are five options: Zstd_BEST_SPEED, Zstd_HIGH_SPEED, Zstd_HIGH_COMPRESSION, Zstd_BEST_COMPRESSION and Zstd_DEFAULT_COMPRESSION. 
                       (Their levels are 1, 3, 19, 22, 3, respectively, default value is Zstd_HIGH_SPEED)

            - randomAccess : supporting Random Access or not. 
                             randomAccess = 1 means that the compression will allow the random access in the decompression
                             Note: need to switch on --enable-randomaccess during the compilation in advance.
                             (Default value is 0)

            [User Parameters <= The following parameters are better to be changed on demand]
             ------------------------------------------------------------------------------
            - errorBoundMode: 8 options to control different types of error bounds (detailed can be found in the user guide of SZ)
                                ABS_AND_REL, ABS_OR_REL, ABS, REL (or VR_REL), PW_REL, ABS_AND_PW_REL, ABS_OR_PW_REL, REL_AND_PW_REL, REL_OR_PW_REL
                                (Default value is REL)

            - absErrBound : absolute Error Bound (NOTE: it's valid when errorBoundMode is related to ABS (i.e., absolute error bound)
                            absErrBound is to limit the (de)compression errors to be within an absolute error. For example, absErrBound=0.0001 means the decompressed value must be in [V-0.0001,V+0.0001], where V is the original true value.
                            (Default value is 4.E-7)

            - relBoundRatio : relative Bound Ratio (NOTE: it's valid only when errorBoundMode is related to REL (i.e., value_range based relative error bound)
                              relErrBound is to limit the (de)compression errors by considering the global data value range size (i.e., taking into account the range size (max_value - min_value)).
                              For example, suppose relBoundRatio is set to 0.01, and the data set is {100,101,102,103,104,...,110}, so the global value range size is 110-100=10, so the error bound will actually be 10*0.01=0.1, from the perspective of "relBoundRatio"

            - psnr : expected PSNR (Note: only valid when errorBoundMode = PSNR)
                     psnr is to spesify the PSNR of the compression. It's valid only when errorBoundMode == PSNR

            - normErr : NORM2 Error: sqrt((x1-x1')^2+(x2-x2')^2+....+(xN-xN')^2) (default value 0.5)

            - pw_relBoundRatio : point-wise relative Bound Ratio (NOTE: only valid when errorBoundMode is related to PW_REL)
                                 pw_relBountRatio is to limit the (de)compression errors by considering the point-wise original data values.
                                 For example, suppose pw_relBoundRatio is set to 0.01, and the data set is {100,101,102,103,104,...,110}, so the compression errors will be limited to {1,1.01,1.02,....1.10} for the data points.
                                 Only valid when errorBoundMode = PW_REL

            - accelerate_pw_rel_compression : superfast compression mode for point-wise relative error bound
                                              True (yes) or False (no)
    )RAW";

PyObject *
py_compress(PyObject *self, PyObject *args)
{
    PyObject *arrays;
    PyObject *options = nullptr;
    if (!PyArg_ParseTuple(args, "O|O", &arrays, &options)) {
        PyErr_SetString(PyExc_SyntaxError,
                        "Wrong syntax. Right syntax : pack(array or list of arrays, dictionnary options or file name");
        return NULL;
    }
    bool                         is_list = false;
    std::vector<PyArrayObject *> np_arrays;
    if (PyList_Check(arrays)) {
        is_list = true;
        np_arrays.reserve(PyList_Size(arrays));
        for (Py_ssize_t i = 0; i < PyList_Size(arrays); ++i)
            np_arrays.push_back((PyArrayObject *)PyList_GetItem(arrays, i));
    }
    else if (PyArray_Check(arrays))
        np_arrays.push_back((PyArrayObject *)arrays);
    else {
        PyErr_SetString(PyExc_TypeError, "First argument must be an array or a list of array");
        return NULL;
    }

    const char *filename = nullptr;
    sz_params parameters;
    init_parameters(parameters);
    if (options) {
        if (PyUnicode_Check(options))
#if PY_VERSION_HEX >= 0x03000000
            filename = PyUnicode_AsUTF8(options);
#else
            filename = PyString_AsString(options);
#endif
        else if (PyDict_Check(options)) {
            bool is_ok = set_parameters_from_dictionnary(options, parameters);
            if (!is_ok) return NULL;
        }
    }
    if (filename) {
        int ok = SZ_Init(filename);
        if (ok == SZ_NSCS) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize SZ librarie ! (Wrong parameters ?)");
            return NULL;
        }
    }
    else {
        int ok = SZ_Init_Params(&parameters);
        if (ok == SZ_NSCS) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize SZ librarie ! (Wrong parameters ?)");
            return NULL;
        }
    }
    PyObject *compressed_list = PyList_New(np_arrays.size());
    int       ind_list        = 0;
    for (const auto &an_array : np_arrays) {
        std::size_t outSize, r5, r4, r3, r2, r1;
        r5              = 0;
        r4              = 0;
        r3              = 0;
        r2              = 0;
        r1              = 0;
        int       ndims = PyArray_NDIM(an_array);
        npy_intp *dims  = PyArray_DIMS(an_array);
        bool is_c_order = false;
        if (PyArray_CHKFLAGS(an_array, NPY_ARRAY_C_CONTIGUOUS)) is_c_order = true;
        if (is_c_order)
        {
            if (ndims > 4) r5 = dims[0];
            if (ndims > 3) r4 = dims[ndims-4];
            if (ndims > 2) r3 = dims[ndims-3];
            if (ndims > 1) r2 = dims[ndims-2];
            r1                = dims[ndims-1];
        }
        else //_ Fortran order _
        {
            r1                = dims[0];
            if (ndims > 1) r2 = dims[1];
            if (ndims > 2) r3 = dims[2];
            if (ndims > 3) r4 = dims[3];
            if (ndims > 4) r5 = dims[4];
        }
        unsigned char *compressed_data = SZ_compress(SZ_DOUBLE, PyArray_DATA(an_array), &outSize, r5, r4, r3, r2, r1);
        npy_intp out_dim = outSize;
        PyObject *shape = PyTuple_New(ndims);
        for (int i = 0; i < ndims; ++i) PyTuple_SET_ITEM(shape, i, PyLong_FromLong(long(dims[i])));
        PyObject *obj = PyTuple_New(3);
        PyTuple_SET_ITEM(obj, 0, shape);
        PyArrayObject* data = (PyArrayObject*)PyArray_SimpleNewFromData(1, &out_dim, NPY_BYTE, compressed_data);
        char* ptr = (char*)PyArray_DATA(data);
        
        PyArray_ENABLEFLAGS(data, NPY_ARRAY_OWNDATA);
        PyTuple_SET_ITEM(obj, 1, (PyObject*)data);
        if (is_c_order)
        {
            Py_IncRef(Py_True);
            PyTuple_SetItem(obj, 2, Py_True); 
        }
        else
        {
            Py_IncRef(Py_False);
            PyTuple_SetItem(obj, 2, Py_False); 
        }
        PyList_SetItem(compressed_list, ind_list, obj);
        ind_list++;
    }
    //display_parameters(parameters);
    SZ_Finalize();
    if (!is_list) {
        PyObject *array = PyList_GetItem(compressed_list, 0);
        Py_INCREF(array);
        Py_DECREF(compressed_list);
        return array;
    }
    return compressed_list;
}

static const char *decompress_doc = R"RAW(
Decompress a compressed array or a list of compressed arrays.
-------------------------------------------------------------
Usage: array = sz.unpack(compressed_array[,options]);
      or
        lst_of_arrays = sz.unpack([comp_arr1, comp_arr2,..., comp_arrb][,options])
      where :
        compressed_array (or comp_arr1 and so.) are compressed data (with sz),
      and options (optional arguments) as same signification than options in compress function.
)RAW";

struct sz_ndims {
    sz_ndims() : r{0, 0, 0, 0, 0} {}
    std::array<std::size_t, 5> r;
    int
    ndims() const
    {
        int n = 5;
        while ((n > 0) && (r[5 - n] == 0)) --n;
        return n;
    }
};

PyObject *
py_decompress(PyObject *self, PyObject *args)
{
    PyObject *cpr_arrays, *options;
    options = nullptr;
    if (!PyArg_ParseTuple(args, "O|O", &cpr_arrays, &options)) {
        PyErr_SetString(PyExc_SyntaxError, "Wrong syntax. Right syntax : unpack(array or list of compressed arrays");
        return NULL;
    }
    //_ _________________ Dépouillement des tableaux à décompresser____________
    bool                         is_list = false;
    std::vector<PyArrayObject *> np_cpr_arrays;
    std::vector<sz_ndims>        shape_arrays;
    std::vector<int>             is_fortran_order;
    if (PyList_Check(cpr_arrays)) {
        is_list              = true;
        Py_ssize_t list_size = PyList_Size(cpr_arrays);
        np_cpr_arrays.reserve(list_size);
        shape_arrays.resize(list_size);
        is_fortran_order.reserve(list_size);
        for (Py_ssize_t i = 0; i < list_size; ++i) {
            PyObject *tuple = PyList_GetItem(cpr_arrays, i);
            if (!PyTuple_Check(tuple)) {
                PyErr_SetString(PyExc_TypeError, "The values of the list must be tuple as (shape, compressed data)");
                return NULL;
            }
            PyObject *array = PyTuple_GetItem(tuple, 1);
            if (!PyArray_Check(array)) {
                PyErr_SetString(PyExc_TypeError, "Second value of the tuple must be an array with compressed data");
                return NULL;
            }
            np_cpr_arrays.push_back((PyArrayObject *)array);
            PyObject *shp = PyTuple_GetItem(tuple, 0);
            if (!PyTuple_Check(shp)) {
                PyErr_SetString(PyExc_TypeError, "A shape must be a tuple of integers");
                return NULL;
            }
            Py_ssize_t dimshape = PyTuple_Size(shp);
            if (dimshape > 5) {
                PyErr_SetString(PyExc_ValueError, "Shape must have only five or less elements");
                return NULL;
            }
            PyObject *py_is_c_order = PyTuple_GetItem(tuple, 2);
            for (Py_ssize_t j = 0; j < dimshape; ++j) {
                PyObject *py_dim = PyTuple_GetItem(shp, j);
                if (PyLong_Check(py_dim) == false && PyInt_Check(py_dim) == false) 
                {
                    PyErr_SetString(PyExc_TypeError, "Values in shape must be integers");
                    return NULL;
                }
                long dim                            = PyLong_AsLong(py_dim);
                if (py_is_c_order == Py_True)
                {
                    shape_arrays[i].r[5 - dimshape + j] = std::size_t(dim);
                }
                else//_ Fortran order _
                {
                    shape_arrays[i].r[5 - j -1] = std::size_t(dim);
                }
            }
            if (py_is_c_order == Py_True)
                is_fortran_order.push_back(0);
            else
            {
                assert(py_is_c_order == Py_False);
                is_fortran_order.push_back(1);
            }
        } // End for (i= 0; i < list_size; ...)
    }     // Fin du cas si c'est une liste
    else if (PyTuple_Check(cpr_arrays)) {
        //_ ___________ Cas où il n'y a qu'un seul tableau ____________________
        PyObject *shape = PyTuple_GetItem(cpr_arrays, 0);
        PyObject *array = PyTuple_GetItem(cpr_arrays, 1);
        if (!PyArray_Check(array)) {
            PyErr_SetString(PyExc_TypeError, "Second value of tuple must be an array with compressed data");
            return NULL;
        }
        np_cpr_arrays.push_back((PyArrayObject *)array);
        shape_arrays.resize(1);
        if (!PyTuple_Check(shape)) {
            PyErr_SetString(PyExc_TypeError,
                            "First value of tuple must be a tuple containing the shape of the original array");
            return NULL;
        }
        Py_ssize_t dimshape = PyTuple_Size(shape);
        if (dimshape > 5) {
            PyErr_SetString(PyExc_ValueError, "Shape must have only five or less elements");
            return NULL;
        }
        PyObject *py_is_c_order = PyTuple_GetItem(cpr_arrays, 2);
        for (Py_ssize_t j = 0; j < dimshape; ++j) {
            PyObject *py_dim = PyTuple_GetItem(shape, j);
            if (PyLong_Check(py_dim) == false && PyInt_Check(py_dim) == false) 
            {
                PyErr_SetString(PyExc_TypeError, "Values in shape must be integers");
                return NULL;
            }
            long dim                            = PyLong_AsLong(py_dim);
            if (py_is_c_order == Py_True)
            {
                shape_arrays[0].r[5 - dimshape + j] = std::size_t(dim);
            }
            else//_ Fortran order _
            {
                shape_arrays[0].r[5 - j -1] = std::size_t(dim);
            }
        }
        is_fortran_order.reserve(1);
        if (py_is_c_order == Py_True)
        {
            is_fortran_order.push_back(0);
        }
        else
        {
            assert(py_is_c_order == Py_False);
            is_fortran_order.push_back(1);
        }
    } // Fin du cas avec un seul tableau
    else {
        PyErr_SetString(PyExc_TypeError, "First argument must be an compressed array or a list of compressed arrays");
        return NULL;
    }

    const char *filename = nullptr;
    sz_params   parameters;
    init_parameters(parameters);
    if (options) {
        if (PyUnicode_Check(options))
#if PY_VERSION_HEX >= 0x03000000
            filename = PyUnicode_AsUTF8(options);
#else
            filename = PyString_AsString(options);
#endif   
        else if (PyDict_Check(options)) {
            bool is_ok = set_parameters_from_dictionnary(options, parameters);
            if (!is_ok) return NULL;
        }
    }
    int ok;
    if (filename) {
        int ok = SZ_Init(filename);
        if (ok == SZ_NSCS) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize SZ librarie ! (Wrong parameters ?)");
            return NULL;
        }
    }
    else {
        int ok = SZ_Init_Params(&parameters);
        if (ok == SZ_NSCS) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize SZ librarie ! (Wrong parameters ?)");
            return NULL;
        }
    }
    //_ ___________Parcours des tableaux pour les décompresser ________________
    PyObject *lst_out_arrays;
    lst_out_arrays = PyList_New(np_cpr_arrays.size());
    for (size_t i = 0; i < np_cpr_arrays.size(); ++i) {
        npy_intp  dims[5];
        int       ndim;
        ndim = shape_arrays[i].ndims();
        bool is_fortran = (is_fortran_order[i]==1);
        if (!is_fortran)
            for (int j = 0; j < ndim; ++j) { dims[j] = shape_arrays[i].r[5 - ndim + j]; }
        else
            for (int j = 0; j < ndim; ++j) { dims[j] = shape_arrays[i].r[4 - j]; }
        for ( int j = ndim; j < 5; ++j) dims[j] = 0;
        PyArrayObject *py_array;
        if (!is_fortran)  py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE,0);
        else 
        {
            py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE,1);
        }
        unsigned char *py_array_data = (unsigned char *)PyArray_DATA(py_array);
        std::size_t    cpr_length    = PyArray_SIZE(np_cpr_arrays[i]);
        unsigned char *cpr_data = (unsigned char *)PyArray_DATA(np_cpr_arrays[i]);
        std::size_t    sz_decompressed_array =
            SZ_decompress_args(SZ_DOUBLE, cpr_data, cpr_length, py_array_data, shape_arrays[i].r[0],
                               shape_arrays[i].r[1], shape_arrays[i].r[2], shape_arrays[i].r[3], shape_arrays[i].r[4]);
        if (!is_list) {
            SZ_Finalize();
            Py_DecRef(lst_out_arrays);
            return (PyObject *)py_array;
        }
        PyList_SetItem(lst_out_arrays, i, (PyObject *)py_array);
    }
    SZ_Finalize();
    return lst_out_arrays;
}
} // namespace sz

static PyMethodDef Pycompressor_sz[] = {{"pack", sz::py_compress, METH_VARARGS, sz::compress_doc},
                                        {"unpack", sz::py_decompress, METH_VARARGS, sz::decompress_doc},
                                        {NULL, NULL}};

#define ADD_INTEGER_KEY(name)                                                                                          \
    {                                                                                                                  \
        PyObject *py_val = PyLong_FromLong(name);                                                                      \
        Py_INCREF(py_val);                                                                                             \
        PyModule_AddObject(m, #name, py_val);                                                                          \
    }


static const char *       module_doc = R"RAW(
Python interface for the scientific dataset compressor sz
---------------------------------------------------------
Two functions to use :
  - pack   : compress an array or a list of arrays containing double values
  - unpack : decompress an array or a list of arrays, providing shapes for decompressed arrays
)RAW";

#if PY_MAJOR_VERSION >= 3
// =====================================================================
static struct PyModuleDef moduledef  = {PyModuleDef_HEAD_INIT,
                                       "csz",
                                       module_doc,
                                       -1, // sizeof(struct module_state),
                                       Pycompressor_sz,
                                       NULL,
                                       NULL, // myextension_traverse,
                                       NULL, // myextension_clear,
                                       NULL};


PyMODINIT_FUNC PyInit_csz(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (m == NULL) return NULL;
    /* Tres important : initialise numpy afin de pouvoir l'utiliser ici !!!! */
    import_array();

    for (const auto &key : sz::sz_keys) {
        PyObject *py_key = PyUnicode_FromString(key);
        Py_INCREF(py_key);
        PyModule_AddObject(m, key, py_key);
    }
    ADD_INTEGER_KEY(ABS);
    ADD_INTEGER_KEY(REL);
    ADD_INTEGER_KEY(ABS_AND_REL);
    ADD_INTEGER_KEY(ABS_OR_REL);
    ADD_INTEGER_KEY(ABS_OR_PW_REL);
    ADD_INTEGER_KEY(ABS_AND_PW_REL);
    ADD_INTEGER_KEY(PW_REL);
    ADD_INTEGER_KEY(PSNR);

    return m;
}
#else
PyMODINIT_FUNC initcsz(void)
{
    PyObject* m = Py_InitModule3("csz", Pycompressor_sz, module_doc);
    if (m == NULL) return;
    /* Tres important : initialise numpy afin de pouvoir l'utiliser ici !!!! */
    import_array();

    for (const auto &key : sz::sz_keys) {
        PyObject *py_key = PyString_FromString(key);
        Py_INCREF(py_key);
        PyModule_AddObject(m, key, (PyObject *)py_key);
    }
}
#endif
