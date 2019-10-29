#ifndef _KCORE_NUMPY_TYPES_HPP_
#define _KCORE_NUMPY_TYPES_HPP_
#include "Numpy/importNumpy.h"
#include <complex>
#include <numpy/arrayobject.h>

namespace KCore {
    namespace Numpy {

        template <typename T>
        struct type {
            static int value() { return NPY_VOID; }
        };

        template <>
        struct type<bool> {
            static int value() { return NPY_BOOL; }
        };

        template <>
        struct type<signed char> {
            static int value() { return NPY_BYTE; }
        };

        template <>
        struct type<unsigned char> {
            static int value() { return NPY_UBYTE; }
        };

        template <>
        struct type<short> {
            static int value() { return NPY_SHORT; }
        };

        template <>
        struct type<unsigned short> {
            static int value() { return NPY_USHORT; }
        };

        template <>
        struct type<int> {
            static int value() { return NPY_INT; }
        };

        template <>
        struct type<unsigned int> {
            static int value() { return NPY_UINT; }
        };

        template <>
        struct type<long> {
            static int value() { return NPY_LONG; }
        };

        template <>
        struct type<unsigned long> {
            static int value() { return NPY_ULONG; }
        };

        template <>
        struct type<long long> {
            static int value() { return NPY_LONGLONG; }
        };

        template <>
        struct type<unsigned long long> {
            static int value() { return NPY_ULONGLONG; }
        };

        template <>
        struct type<float> {
            static int value() { return NPY_FLOAT32; }
        };

        template <>
        struct type<double> {
            static int value() { return NPY_DOUBLE; }
        };

        template <>
        struct type<std::complex<float>> {
            static int value() { return NPY_COMPLEX64; }
        };

        template <>
        struct type<std::complex<double>> {
            static int value() { return NPY_COMPLEX128; }
        };

        template <>
        struct type<PyObject> {
            static int value() { return NPY_OBJECT; }
        };
    } // namespace Numpy
} // namespace KCORE

#endif
